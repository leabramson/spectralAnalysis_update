pro pyspec1ID, ID, PA, z

  ID = strcompress(string(ID), /rem)
  PA = strcompress(string(PA), /rem)
  z = strcompress(string(z), /rem)
  
  print, ''
  print, 'Running PYSPECFIT for ID '+ID+' PA '+Pa
  print, ''

  file_mkdir, ID+'_results'
  spawn, 'mpirun -np 12 python fitspec.py '+ID+'_'+PA+' OPTFL '+z+' 0.05 '+ID+'_results/optfl'+PA
  spawn, 'mpirun -np 12 python fitspec.py '+ID+'_'+PA+' INNFL '+z+' 0.05 '+ID+'_results/innfl'+PA
  spawn, 'mpirun -np 12 python fitspec.py '+ID+'_'+PA+' INTUP '+z+' 0.05 '+ID+'_results/intup'+PA
  spawn, 'mpirun -np 12 python fitspec.py '+ID+'_'+PA+' INTDN '+z+' 0.05 '+ID+'_results/intdn'+PA
  spawn, 'mpirun -np 12 python fitspec.py '+ID+'_'+PA+' OUTUP '+z+' 0.05 '+ID+'_results/outup'+PA
  spawn, 'mpirun -np 12 python fitspec.py '+ID+'_'+PA+' OUTDN '+z+' 0.05 '+ID+'_results/outdn'+PA

end
