import os

cp=os.getcwd()
for i in os.listdir():
    target_p=os.path.join(cp,i)
    try:
        os.chdir(target_p)
        os.chdir('relax')
        os.system('sbatch submit.sh')
        os.chdir('..')
        os.chdir('relax_xtal')
        os.system('sbatch submit.sh')
        os.chdir(cp)
    except:
        print('This target fail')
        os.chdir(cp)
    
    


    