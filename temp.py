from simuojo import Structure

# CTCTGTTAGTATTTTCTGCGTGCCAGTGAGTCGGATCTCCACAGTTCTCTGTG
a=Structure('GCATUGCATU',name='new')
a.fold()

_=a.plot_fold(save=False)

s0=a.dotgraph[0]
s0.plot_2d()

a._dot

a._pairtuple
