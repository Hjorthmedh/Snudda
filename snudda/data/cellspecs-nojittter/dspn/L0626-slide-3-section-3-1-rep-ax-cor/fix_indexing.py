
#swcFile = 'WT-0728MSN01-cor-rep-ax-res5.swc'
swcFile = 'WT-0728MSN01-cor-rep-ax.swc'
with open(swcFile,'r') as f:
    lines = f.readlines()

to_write = ''
c = 1
a0 = 0
mapper = {}
for ss in lines:
    if(ss[0] != '#'):
        values = [s for s in ss.split()]
        if not int(values[1]) == 2:
            to_write = '{}{}'.format(to_write, ss)
        else:
            if not a0:
                a0 = int(values[0]) - c
            else:
                if values[6] == aprev: values[6] = str( c-1 )
                else: 
                    print( values[6],values[0] )
                    print( mapper[values[6]]-1 )
                    values[6] = str( mapper[values[6]] )
            aprev = values[0]
            mapper[values[0]] = c
            values[0] = str(c)
            to_write = '{}{}\n'.format(to_write, ' '.join(values))
            
        c += 1  
    
with open('t.swc','w') as f:
    f.write(to_write)    
