TITLE AMPACOD  

COMMENT
	AMPA model from Paper (version 15 sept 2004).
	
ENDCOMMENT

NEURON {
	POINT_PROCESS AmpaCOD	
	NONSPECIFIC_CURRENT i	
	RANGE r1FIX,r2,r3,r4,r5,r1,r6,r6FIX
	RANGE g,gmax,kB,Cdur,Erev 
	RANGE gg1,gg2,gg3,Tdiff
	RANGE T,Tmax,Trelease 	
	RANGE tau_1,tau_rec,tau_facil,U,u0	
	RANGE A1,A2,A3,tau_dec1,tau_dec2,tau_dec3		: comes from fit 
	RANGE tdelay,ton	 
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)	
	(mM) = (milli/liter)
	(pS) = (picosiemens)
	(nS) = (nanosiemens)
	(um) = (micrometer)
	PI   = (pi)(1)
}

PARAMETER {
	: Parametri Postsinaptici
	r1FIX		= 5.4		(/ms/mM) 	 
	r2		= 0.82		(/ms)		 
	r3		= 0		(/ms)		 
	r4		= 0		(/ms)		 
	r5		= 0.013		(/ms)		
	r6FIX		= 1.12		(/ms/mM)	
	gmax		= 700 		(nS)		 
	Cdur		= 0.3		(ms)		 
	Erev		= 0		(mV)
	kB		= 0.44		(mM)		
	: Diffusion: M=21500, R=1.033, D=0.223, lamd=0.02			
	A1 			= 0.131837 
	A2			= 0.0555027	 
	A3 			= 0.0135232	 
	tau_dec1 		= 3.4958	 
	tau_dec2 		= 16.6317	 
	tau_dec3 		= 128.983	

	: Parametri Presinaptici
	tau_1 		= 3 (ms) 	< 1e-9, 1e9 >
	tau_rec 	= 35.1 (ms) 	< 1e-9, 1e9 > 	
	tau_facil 	= 10.8 (ms) 	< 0, 1e9 > 	

	U 		= 0.416 (1) 	< 0, 1 >
	u0 		= 0 (1) 	< 0, 1 >	: se u0=0 al primo colpo y=U
	Tmax		= 1  (mM)
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: conductance
	r1		(/ms)
	r6		(/ms)
	T		(mM)
	Trelease	(mM)
	Tdiff		(mM)
	tdelay		(ms)
	Tdiff_0		(mM)
	ton		(ms)	
	
	x
	PRE
}

STATE {	
	C
	O
	D
	gg1
	gg2
	gg3
	sink
}	
	

INITIAL {
	C		=	1
	O		=	0
	D		=	0
	T		=	0 	(mM)
	Tdiff		=	0	(mM)
	Trelease	=	0 	(mM)
	Tdiff_0		=	0	(mM)
	gg1		=	0
	gg2		=	0
	gg3		=	0   
	ton		=  	-1   (ms)
	PRE		=	0
}
FUNCTION SET_tdelay(R,D){ tdelay=0.25*R*R/D } 

BREAKPOINT {
	if( (t-ton)>tdelay  ){
		Tdiff=gg1+gg2+gg3
		Tdiff_0 = Tdiff
	}else{
		Tdiff=Tdiff_0+(A1+A2+A3)*PRE*(t-ton)/tdelay
	}
	Trelease=T+Tdiff
	SOLVE kstates METHOD sparse
	g =gmax * O
	i = (1e-6) * g * (v - Erev) 
}


KINETIC kstates {
	: Postsynaptic scheme
	r1 = r1FIX * Trelease^2 / (Trelease + kB)^2
	r6 = r6FIX * Trelease^2 / (Trelease + kB)^2
	~ C  <-> O	(r1,r2)
	~ O  <-> D	(r3,r4)
	~ D  <-> C	(r5,r6)
	CONSERVE C+O+D = 1
	: Glutamate diffusion wave
	~ gg1 <-> sink (1/tau_dec1,0)
	~ gg2 <-> sink (1/tau_dec2,0)
	~ gg3 <-> sink (1/tau_dec3,0)
}


NET_RECEIVE(weight, on, nspike, flagtdel,t0 (ms),y, z, u, tsyn (ms)) {
	INITIAL {
		flagtdel=1
		nspike = 1
		Tdiff=0
		y=0
		z=0
		u=u0
		tsyn=t
	}
   	if (flag == 0) { 
		nspike = nspike + 1
		if (!on) {
			ton=t
			t0=t
			on=1				
			z=z*exp(-(t-tsyn)/tau_rec)
			z=z+(y*(exp(-(t - tsyn)/tau_1)-exp(-(t-tsyn)/tau_rec))/(tau_1/tau_rec-1))
			y=y*exp(-(t-tsyn)/tau_1)			
			x=1-y-z
			if(tau_facil>0){ 
				u=u*exp(-(t-tsyn)/tau_facil)
				u=u+U*(1-u)							
			}else{u=U}
			y=y+x*u
			PRE=y
			T=Tmax*y
			tsyn=t						
		}
		net_send(Cdur,nspike)
		net_send(tdelay,flagtdel)						
    	}
	if(flag==nspike){ 		
		T = 0
		on = 0
	}
	if (flag == flagtdel){
		flagtdel = flagtdel+1
		state_discontinuity(gg1,gg1+A1*x*u)	 
		state_discontinuity(gg2,gg2+A2*x*u)	 
		state_discontinuity(gg3,gg3+A3*x*u)	 
	}
}	 

 
