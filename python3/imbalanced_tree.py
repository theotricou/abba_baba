# Theo
# create completly imbalance phylo trees

from ete3 import Tree as tr

t1=tr()
t1.populate(2, random_branches=False)

for n in range(1,99):
    node=None
    for i in t1:
        if i.name == "aaaaaaaaaa":
            node=i
            t2=tr()
            t2.populate(2, random_branches=False)
            t2.name="t2"
            break
    node.add_child(t2, dist=0)
    rem=t1&"t2"
    rem.delete()
    t1.convert_to_ultrametric()
    t=t1.copy()
    names=1
    for i in t:
        i.name=names
        names+=1
    t.write(format=1, outfile=str(n+2), format_root_node=True)
