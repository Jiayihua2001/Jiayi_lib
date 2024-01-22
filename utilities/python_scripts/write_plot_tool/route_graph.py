from graphviz import Digraph

dot = Digraph(comment='Dimer Generation Process')
dot.attr(size='6,6', ratio='fill', bgcolor='lightgrey')
dot.node_attr.update(style='filled', color='black', fillcolor='skyblue', shape='box', fontname='Helvetica', fontsize='12')
dot.edge_attr.update(color='blue', arrowhead='vee', fontsize='10')


dot.node('A', 'Input mol1 mol2')
dot.node('B', 'Central Molecule Selection')
dot.node('C', 'Grid generation')
dot.node('D', 'Dimer generation and validation')
dot.node('E', 'Duplication removal')
dot.node('F', 'Optimization(offset)')
dot.node('G', 'Select Stable Oritation')
dot.node('H', 'Relaxation')

dot.edge('A', 'B', label='Depend on mol symmetry')
dot.edge('B', 'C', label='Box Size Detection and Adjustment')
dot.edge('C', 'D', label='Distance calculation and check')
dot.edge('D', 'E', label='')
dot.edge('E', 'G', label='Calculate the energy')
dot.edge('G', 'H', label='for instance: select top 100 most stable structures for relaxation')
dot.edge('E', 'F', label='No change on orientation,only distance')


dot.render('process_graph.gv', view=True)
