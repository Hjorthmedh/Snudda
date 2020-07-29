
swcFile = 'Experiment-9862corrected-cor-rep.swc'
with open(swcFile,'r') as f:
    lines = f.readlines()

to_write = ''
c = 1

mapper = {}
for ss in lines:
    if(ss[0] != '#'):
        values = [s for s in ss.split()]
        if values[6] == '-1':
            org = values[2:5]
            for i in range(2,5): values[i] = '0'
        else:
            for i in range(2,5):
                values[i] = '{:.3f}'.format( float(values[i])-float(org[i-2]) )
        to_write = '{}{}\n'.format(to_write, ' '.join(values))
        
with open('t.swc','w') as f:
    f.write(to_write)    
