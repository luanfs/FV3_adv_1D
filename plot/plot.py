import numpy as np
import matplotlib.pyplot as plt

# Python script to plot the outputs

# directory
graphdir ='../graphs/' # must exist
datadir  ='../data/' # must exist
figformat = 'png'

# some constants
N    = 48 # number of cells
tc   = 2  # test case
hord = 8  # advection scheme
nplots = 7

# basename for plotting
basename = "tc"+str(tc)+"_N"+str(N)+"_hord"+str(hord)+"_t"

# x axis points for plotting
x = np.linspace(0, 1, N)

for t in range(0, nplots+1):
   input_name  = datadir+basename+str(t)+'.txt'
   output_name = graphdir+basename+str(t)+'.'+figformat

   data = np.loadtxt(input_name)

   time = data[0]
   time = str("{:.2e}".format(time))

   massvar = data[1]
   massvar = str("{:.2e}".format(massvar))

   cfl = data[2]
   cfl = str("{:.2e}".format(cfl))

   y = data[3:]

   # plot the graph
   plt.plot(x,y)

   if tc == 1:
       plt.ylim(-0.1, 1.1)
   elif tc == 2:
      plt.ylim(-0.1, 2)
 
   # Label
   plt.xlabel('$x$')
   plt.ylabel('$y$')
   #plt.legend()
   plt.grid(True, which="both")
   title = "N="+str(N)+", hord="+str(hord)+", time = "+time+" days"+\
   "\nCFL="+cfl+", mass variation="+massvar  
   plt.title(title)
   plt.savefig(output_name+'.'+figformat, format=figformat)
   plt.close()
   print("Plotted "+ output_name)
   

