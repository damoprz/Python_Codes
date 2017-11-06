import numpy as np
import read_muram as rmu
import dp_plot_tools as dplt

RUNDIR = '/home/przybylski/scratch/3D_TEST_HIGH/'

snap = rmu.MURaM_snap(RUNDIR,boxtop=12.0e8)
snap.load('190000',tooload=['rho','tem','vz','bz','Qtot','vx','vy','bx','by','tau','eps','pre']) ##,'ne','QxH','QxCa','QxMg','QxCor'])

scale = (snap.tau*snap.tau/(snap.tau*snap.tau+1.0e-16))

tau_av = np.mean(snap.tau,axis=(0,1))
for i in range(snap.nz):
   if (tau_av[i] <= 1):
     z_tau_cor = snap.zax[i-1] + (snap.zax[i]-snap.zax[i-1])/(tau_av[i] - tau_av[i-1])*(1.0-tau_av[i-1])
     print(snap.zax[i],tau_av[i],z_tau_cor)
     snap.zax = snap.zax-z_tau_cor
     break

xr = [snap.gxmin[0]/1.0e8,snap.gxmax[0]/1.0e8]
yr = [snap.gxmin[1]/1.0e8,snap.gxmax[1]/1.0e8]
zr = [snap.gxmin[2]/1.0e8,snap.gxmax[2]/1.0e8]-z_tau_cor/1.0e8

xyslice = 228
xzslice = 95

z_plotr = [-0.0,3.0]
hr = [1.15,1.0]
asp = [1,'auto']
fontsize = 12

[fig_vh,gs_vh] = dplt.fig_open(figsize=[17.5,6.8],numx=4,numy=2,hr=hr)

[im1,im2,cbt1,cbt2]=dplt.plotvhslice(fig_vh,snap.vz[:,:,:]*1.0e-5,xyslice,xzslice,gs_vh[0],gs_vh[4],1.0e-8*snap.xax,1.0e-8*snap.yax,1.0e-8*snap.zax,xr,yr,z_plotr,br=[-10.0,10.0],cmap='seismic_r',cbar='Right',sym=True,xtitle=True,ytitle=True,title=r'$v_z (km\; s^{-1})$',asp=asp,fontsize=fontsize)

[im1,im2,cbt1,cbt2]=dplt.plotvhslice(fig_vh,snap.bz[:,:,:],xyslice,xzslice,gs_vh[1],gs_vh[5],1.0e-8*snap.xax,1.0e-8*snap.yax,1.0e-8*snap.zax,xr,yr,z_plotr,br = [-50,50,-200,200],cmap='seismic_r',cbar='Right',sym=True,xtitle=True,ytitle=False,title=r'$B_z$ (Gauss)',asp=asp,fontsize=fontsize)

[im1,im2,cbt1,cbt2]=dplt.plotvhslice(fig_vh,np.log10(snap.rho[:,:,:]),xyslice,xzslice,gs_vh[2],gs_vh[6],1.0e-8*snap.xax,1.0e-8*snap.yax,1.0e-8*snap.zax,xr,yr,z_plotr,cbar='Right',sym=False,xtitle=True,ytitle=False,title=r'$\log10(\rho)\; (g\;cm^{-3})$',asp=asp,fontsize=fontsize)

[im1,im2,cbt1,cbt2]=dplt.plotvhslice(fig_vh,snap.tem,xyslice,xzslice,gs_vh[3],gs_vh[7],1.0e-8*snap.xax,1.0e-8*snap.yax,1.0e-8*snap.zax,xr,yr,z_plotr,br = [3.0e3,1.0e4],cmap='YlOrRd_r',cbar='Right',sym=False,xtitle=True,ytitle=False,title=r'$T$ (K)',asp=asp,fontsize=fontsize)

gs_vh.tight_layout(fig_vh,h_pad=0.5,w_pad=0.5,pad=0.3, rect=[0.0,0.0,1.0,1.0])
fig_vh.show()

