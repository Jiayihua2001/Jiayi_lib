import logging
import os
import sys
import time

from mpi4py import MPI

from gnrs.parser import UserSettingsParser
from gnrs.parser import UserSettingsSanityChecker
from gnrs.core.logging import GenarrisLogger
from gnrs.core import folders
from gnrs.parallel.test import test_bcast
from gnrs.parallel.structs import DistributedStructs
from gnrs.parallel import init_parallel
import gnrs.restart as restart
import gnrs.output as gout
from ase.io import read
# uncomment for  dimer csp
import json


# For global exception handle
comm = MPI.COMM_WORLD


def run():
    """
    When calling genarris_master.py,
    should specifiy the path to an instruction file
    """
    return Genarris(sys.argv)


class Genarris:
    """
    Defines the flow of control in Genarris
    """

    def __init__(self, argv):
        """Defines abstract workflow of Genarris."""

        # The two important member variables
        self.user_settings = dict()
        self.runtime_settings = dict()
        self.restart = False

        self._mpi_init()
        self._log_init()
        self._output_init()
        self._parallel_init()
        self._runtime_settings_init(argv)
        self._restart_init()
        self.get_user_settings(argv)  # Parse input into user_settings
        self.attempt_restart()
        self._folders_init()
        self.comm.barrier()  # Wait for folder creation
        self.run_genarris()
        return

    def _log_init(self):
        """Initialize logger; pass communicator to define master for log"""
        self.Genlogger = GenarrisLogger(self.comm)
        self.logger = logging.getLogger("master")

    def _mpi_init(self):
        """Initialise MPI"""
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.is_master = True if self.rank == 0 else False

    def _parallel_init(self):
        init_parallel(self.comm)

    def _restart_init(self):
        restart.restart_init(self.comm, self.user_settings, self.runtime_settings)

    def _output_init(self):
        """Print hostname, citation, version, parallelization"""
        gout.init_output(self.comm)
        gout.welcome_message()

    def get_user_settings(self, argv):
        """
        Gets the input settings from the config file.

        If gnrs is launched with 'restart' argument instead of a config,
        user_settings from the restart.json are used.

        Otherwise, config file is parsed to build the user_settings
        dict.
        """
        # Hard restart
        if argv[-1] == "restart":
            self.restart = True
            gout.emit("Skip parsing config file. Forcing restart...")
            return

        restart_flag = None
        self.config_path = os.path.abspath(argv[-1])
        self.runtime_settings["config_path"] = self.config_path
        
        # Parse Config file
        gout.print_title("Parsing User Input")
        gout.emit(f"Reading {self.config_path}.")
        # parse config file
        if self.is_master:
            parser = UserSettingsParser(self.config_path)
            user_settings = parser.load_user_settings()
            # Update log level
            new_level = user_settings["master"]["log_level"]
            self.Genlogger.reset_loglevel(new_level)
            UserSettingsSanityChecker(user_settings)
            self.user_settings.update(user_settings)
            restart_flag = user_settings["master"]["restart"]

        # Check if restart is requested in ui.conf. If so, set
        # self.restart = True to attempt restart.
        restart_flag = self.comm.bcast(restart_flag, root=0)
        if restart_flag:
            gout.emit("Restart requested in config file." " Attempting restart...")
            if restart.check_restart():
                self.restart = True
                gout.emit("Restart Success!")
            else:
                gout.emit("Restart Failed." " Starting Genarris using config file.")
                self.restart = False

        # Send the config data to other processes
        self.user_settings = self.comm.bcast(self.user_settings, root=0)
        return
    
    ## Add some extra settings for dimer csp


    def dimer_setting(self,user_settings,num_1,num_2):

        out_dir=user_settings
        mol_path_l=out_dir["rtm_set"]["molecule_path"]
        geo_p=mol_path_l[0]
        head, tail = os.path.split(geo_p)
        mol1_p=os.path.join(head,'mol1.in')
        mol2_p=os.path.join(head,'mol2.in')
        mol_path_l.remove(geo_p)
        mol_path_l.append(mol1_p)
        mol_path_l.append(mol2_p)
        out_dir["rtm_set"]["molecule_path"]=mol_path_l
        n_mol_l=out_dir["rtm_set"]["n_atoms_in_mol"]
        n_dimer=n_mol_l[0]
        try:
            if num_1 + num_2 == n_dimer:
                n_mol_l.remove(n_dimer)
                n_mol_l.append(num_1)
                n_mol_l.append(num_2)
        except:
            print(f'error: total n of atoms {n_dimer} != that sum of mol1 and mol2: {num_1}+{num_2}')
        out_dir["rtm_set"]["n_atoms_in_mol"]=n_mol_l
        
        out_dir["usr_set"]["generation"]["generation_type"]='cocrystal'
        out_dir["usr_set"]["generation"]["stoichiometry"]= [1,1]
        out_dir["usr_set"]["master"]["molecule_path"] = ["./mol1.in","./mol2.in"]

        user_setting= out_dir

        return user_setting

    def _runtime_settings_init(self, argv):
        """
        Add to runtime_settings dict:
            config path
            genarris start time
            number of ranks
        Following dictionaries are also attached
            Completion status of tasks
            task start time
            task end time
        """
        self.logger.info("Setting runtime values")
        # Set work dir, temp dir and struct dir
        self.runtime_settings["work_dir"] = self.work_dir = os.getcwd()
        struct_dir = os.path.join(self.work_dir, "structures")
        tmp_dir = os.path.join(self.work_dir, "tmp")
        self.runtime_settings["struct_dir"] = struct_dir
        self.runtime_settings["tmp_dir"] = tmp_dir
        self.runtime_settings["energy_list"] = []
        self.runtime_settings["genarris_start_time"] = time.time()
        self.runtime_settings["size"] = self.size
        return

    def attempt_restart(self):
        if not self.restart:
            return
        restart.load_restart()
        gout.print_title("restarting genarris")
        gout.print_user_settings(self.user_settings)
        gout.double_separator()

    def _folders_init(self):
        """
        Creates tmp and structures directory in the working directory
            if it doesn't exist.
        Creates a molecule folder in tmp and copies molecule into it.
        """
        folders.init_folders(self.is_master)
        if not self.restart:
            self.logger.info("Setting up folders: structures and tmp")
            folders.setup_main_folders(self.runtime_settings)
            folders.copy_molecule(self.user_settings, self.runtime_settings)
        return

    def run_genarris(self):

        from gnrs.generation import StructureGenerationTask
        from gnrs.descriptor import DescriptorEvaluationTask
        from gnrs.cluster import ClusterSelectionTask
        from gnrs.energy import EnergyCalculationTask
        from gnrs.optimize import GeometryOptimizationTask

        self.logger.info("Starting genarris tasks")

        if not restart.is_task_completed("generation"):
            sg = StructureGenerationTask(
                self.comm, self.user_settings, self.runtime_settings
            )
            sg.run()
            restart.write_restart()

        else:
            self.logger.info("Generation task was completed before restart")
            gout.skip_task("generation")

        test_bcast()
        self._check_if_exp_found()


        if self.is_master:
            os.system('cp ./mol1.in ./tmp/molecule')
            os.system('cp ./mol2.in ./tmp/molecule')
            mol_1=read('./mol1.in')
            mol_2=read('./mol2.in')
            os.system('rm -f ./tmp/molecule/geometry_0.in')
            restart_path = os.path.join('./tmp', "restart.json")
            parser = UserSettingsParser(restart_path)
            user_settings = parser.load_user_settings()
            num_1=len(mol_1)
            num_2=len(mol_2)
            user_settings=self.dimer_setting(user_settings,num_1,num_2)
            # Update log level
            self.logger.info("Change user settings for dimer rigid_press")
            UserSettingsSanityChecker(user_settings)
            self.user_settings=user_settings['usr_set']
            self.runtime_settings = user_settings['rtm_set']
        self.user_settings = self.comm.bcast(self.user_settings, root=0)
        self.runtime_settings = self.comm.bcast(self.runtime_settings, root=0)

        self._restart_init()

        if not restart.is_task_completed("rigid_press"):
            g_opt = GeometryOptimizationTask(
                self.comm, self.user_settings, self.runtime_settings, "rigid_press"
            )
            g_opt.run()
            restart.write_restart()
        else:
            self.logger.info("Rigid press optimization completed before restart")
            gout.skip_task("rigid_press")

        test_bcast()
        self._check_if_exp_found()

        # if not restart.is_task_completed("dftbp"):
        #     e_calc = EnergyCalculationTask(
        #         self.comm, self.user_settings, self.runtime_settings, "DFTBP"
        #     )
        #     e_calc.run()
        #     restart.write_restart()
        # else:
        #     gout.skip_task("DFTBP")

        # test_bcast()

        # if not restart.is_task_completed("nocluster-min"):
        #     noclstr = ClusterSelectionTask(
        #         self.comm, self.user_settings, self.runtime_settings, "nocluster"
        #     )
        #     noclstr.run()
        #     restart.write_restart()
        # else:
        #     gout.skip_task("nocluster")

        # if not restart.is_task_completed("acsf"):
        #     det = DescriptorEvaluationTask(
        #         self.comm, self.user_settings, self.runtime_settings, "ACSF"
        #     )
        #     det.run()
        #     restart.write_restart()

        # else:
        #     self.logger.info("ACSF task was completed before restart")
        #     gout.skip_task("acsf")

        # if not restart.is_task_completed("kmeans-center"):
        #     cst = ClusterSelectionTask(
        #         self.comm, self.user_settings, self.runtime_settings, "kmeans"
        #     )
        #     cst.run()
        #     restart.write_restart()
        # else:
        #     self.logger.info("kmeans selection task was completed before restart")
        #     gout.skip_task("kmeans-center")

        # test_bcast()

        # if not restart.is_task_completed("lbfgs-dftbp"):
        #     g_opt = GeometryOptimizationTask(
        #         self.comm, self.user_settings, self.runtime_settings, "lbfgs", "dftbp"
        #     )
        #     g_opt.run()
        #     restart.write_restart()
        # else:
        #     self.logger.info("lbfgs-dftbp completed before restart")
        #     gout.skip_task("lbfgs-dftbp")
        # """
        # if not restart.is_task_completed("dftbp_lbfgs"):
        #     g_opt = GeometryOptimizationTask(
        #         self.comm, self.user_settings, self.runtime_settings, "dftbp_lbfgs"
        #     )
        #     g_opt.run()
        #     restart.write_restart()
        # """
        # self._check_if_exp_found()

        return
        raise SystemExit()

    def _check_if_exp_found(self):
        from ase.io import read
        from gnrs.parallel.io import read_parallel

        try:
            exp = read(self.user_settings["experimental_structure"]["path"], parallel=False)
        except KeyError:
            return
        gout.emit("Searching for Experimental structure within the pool...")
        spath = self.runtime_settings["last_struct_path"]
        structs = read_parallel(spath)
        match_list = DistributedStructs(structs).perform_xtal_match(exp)
        if match_list:
            gout.emit("Found Experimental structure within the pool.")
        else:
            gout.emit("Experimental structure not found within the pool.")
        return


def main():
    try:
        run()
    except:
        logger = logging.getLogger("main")
        logger.exception("Genarris exiting due to error")
        gout.emit(
            "\n***Genarris has abruptly exited." " Please see the log for Traceback***"
        )
        comm.Abort()
    else:
        gout.emit("All tasks completed.\nHave a nice day! :-)")
    return


if __name__ == "__main__":
    main()
