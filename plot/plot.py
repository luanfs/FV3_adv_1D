import numpy as np
import matplotlib.pyplot as plt

# Python script to plot the outputs

# directory
graphdir ='../graphs/' # must exist
datadir  ='../data/' # must exist
figformat = 'png'

# some constants
N    = 48	 # number of cells
tc   = 1  # test case
hords = (8,0)  # advection schemes
dps   = (1,1)  # departure point schemes
nplots = 13

# domain length
erad = 6371.0 # earth radius (km)
Lx = 0.5*np.pi*erad
Lxo2 = Lx*0.5
 
# x axis points for plotting
x = np.linspace(-Lxo2, Lxo2, N)

for t in range(0, nplots+1):
   if t>0:
       # plot ic
       plt.plot(x, y0, color = 'black', label = 'IC')
   for m in range(0, len(dps)):
      dp = dps[m]
      hord = hords[m]

      # basename for plotting
      basename = "tc"+str(tc)+"_N"+str(N)+"_hord"+str(hord)+"_dp"+str(dp)+"_t"
      input_name  = datadir+basename+str(t)+'.txt'
      output_name = graphdir+'adv1d_'+basename+str(t)+'.'+figformat

      data = np.loadtxt(input_name)

      # plot the graph
      y = data[3:]
      plt.plot(x, y, label = 'hord'+str(hord)+'.dp'+str(dp))

   time = data[0]
   time = str("{:.2e}".format(time))

   massvar = data[1]
   massvar = str("{:.2e}".format(massvar))

   cfl = data[2]
   cfl = str("{:.2e}".format(cfl))


   if tc == 1 or tc == 2:
      plt.ylim(-0.2, 1.2)
   elif tc == 3:
      plt.ylim(-0.1, 2.2)
 
   if t == 0:
      y0 = y

   # Label
   plt.xlabel('$x$ (km)')
   plt.ylabel('$y$')
   plt.legend()
   plt.grid(True, which="both")
   title = "N="+str(N)+", time = "+time+" days, CFL="+cfl#+", mass variation="+massvar  
   plt.title(title)
   plt.savefig(output_name, format=figformat)
   plt.close()
   print("Plotted "+ output_name)
   

