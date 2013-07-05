TITLE 

COMMENT
	Reference: Pugh JR and Raman I, Biophysical Journal Volume 88, March 2005 1740-1754
	Model adapted from patch to slice.
ENDCOMMENT

NEURON {
	POINT_PROCESS GRC_GABA
	NONSPECIFIC_CURRENT i
	RANGE g,Cdur,Erev,Open,OpenScaled,ScaleFactor

	:RANGE r1,r2,kon,koff,d1,d2,a1,a2,b1,b2
	RANGE kon,koff,d3,r3,d1d2,r1r2,a1,b1,a2,b2,r1,r2,d1,d2

	RANGE Tmax,gmax,onSET	
	
	RANGE tau_1,tau_rec,tau_facil,U,T 
	RANGE diff_flag,M,Rd,Diff,lamd
	RANGE diffusione,nd	
}

UNITS {
	(nA) 	= (nanoamp)
	(mV) 	= (millivolt)
	(umho)  = (micromho)
	(mM) 	= (milli/liter)
	(pS) 	= (picosiemens)
	PI   	= (pi)(1)
}

PARAMETER {
	: Parametri Postsinaptici
	gmax	=  756.35	(pS)	:1750 
	Cdur	= 0.3		(ms)	

	kon	= 20		(/ms/mM)	
	koff	= 2		(/ms) 		 
	d3	= 15		(/ms) 		 
	r3	= 3.75		(/ms) 		: 0.15, use 3.75 for slices 

	d1d2	= 15		(/ms/mM)	
	r1r2	= 0.007		(/ms)

	a1	= 0.06		(/ms)
	b1	= 0.03		(/ms)
	a2	= 0.4		(/ms)
	b2	= 10		(/ms)
	
	r1	= 7e-4		(/ms)
	r2	= 6e-3		(/ms)
	d1	= 3.3e-4	(/ms)
	d2	= 1.2		(/ms)

	Erev	= -65	(mV)

	: Parametri Presinaptici
	tau_1 		= 0.1 (ms) 	< 1e-9, 1e9 >
	tau_rec 	= 43.4 (ms) 	< 1e-9, 1e9 > 	:55.11 15.7	(first fit!)
	tau_facil 	= 6.22 (ms) 	< 0, 1e9 >    	:2.66 4.85  (first fit!)	
	U 		= 0.35		< 0, 1 >	:0.24  0.18	(first fit!)
	
	Tmax	= 1  (mM)	
	onSET	= 1
		
	: Diffusion parameters
	: Diffusion: M=21.500, R=1.033, D=0.223, lamd=0.02 as in excitatory synapses	

	M		= 52.76		:	46.93			: 20.95 (first fit!)
: numero di (kilo) molecole in una vescicola		
	Rd		= 4.79 (um)	:4.96	: 4.96 (first fit!)
	Diff		= 0.223 (um2/ms)
	lamd		= 20 	(nm)
	diff_flag	= 1			: flag diffusion on/off
	nd		= 1			: kernel exponent of diffusion

	ScaleFactor	= 1 			: for fit purposes
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: conductance
	Open
	OpenScaled	: for fit purposes
	
	T		(mM)	
	Trelease	(mM)
	Mres		(mM)	
	tpre		(ms)

	tspike[50]	(ms)	: will be initialized by the pointprocess
	PRE[50]
	numpulses
	tzero
}

STATE {	
	C
	CA1
	CA2
	DA1
	DA2
	DA2f
	OA1
	OA2	
}

INITIAL {
	C=1
	CA1=0
	CA2=0
	DA1=0
	DA2=0
	DA2f=0
	OA1=0  	
	OA2=0
	CA1=0
	CA2=0
	Open=0
	T=0 		(mM)
	:tpre=1e8	(ms)
		
	numpulses=0
	Mres=1e3* (1e3 * 1e15 / 6.022e23 * M)     : (M) to (mM) so 1e3, 1um^3=1dm^3*1e-15 so 1e15
	FROM i=1 TO 50{ PRE[i-1]=0 tspike[i-1]=0}
	tspike[0]=1e12	(ms)
	if(tau_1>=tau_rec){ 
		printf("Warning: tau_1 (%g) should never be higher neither equal to tau_rec (%g)!\n",tau_1,tau_rec)
		tau_rec=tau_1+1e-5
		:printf("tau_rec has been set to %g\n",tau_rec) 
	} 
}

FUNCTION diffusione(){
	LOCAL DifWave,i	
	DifWave=0
	FROM i=1 TO numpulses{
		tzero=tspike[i-1]
		if(t>tzero){
			DifWave=DifWave+PRE[i-1]*Mres*exp(-Rd*Rd/(4*Diff*(t-tzero)))/((4*PI*Diff*(1e-3)*lamd)*(t-tzero))^nd
		}
	}	
	diffusione=DifWave :Mres*exp(-Rd*Rd/(4*Diff*(t-tpre)))/((4*PI*Diff*(1e-3)*lamd)*(t-tpre))	
}


BREAKPOINT {
	SOLVE kstates METHOD sparse
	Open = OA1 + OA2
	OpenScaled=Open*ScaleFactor
	g = gmax * Open
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {
	if ( diff_flag ) { Trelease = T + diff_flag * diffusione() } else { Trelease = T }
	: second row
	~	C  	<-> 	CA1	(2*kon*Trelease,koff)
	~	CA1 	<-> 	CA2	(kon*Trelease,2*koff)
	~	CA2	<->	DA2f	(d3,r3)
	: third row
	~ 	DA1  	<-> 	DA2	(d1d2*Trelease,r1r2)
	: first <=> second row
	~ 	OA1  	<-> 	CA1	(a1,b1)
	~ 	OA2  	<-> 	CA2	(a2,b2)
	: third <=> second row
	~	DA1	<->	CA1	(r1,d1)
	~	DA2	<->	CA2	(r2,d2)
	CONSERVE C+CA1+CA2+DA1+DA2+DA2f+OA1+OA2 = 1
}


NET_RECEIVE(weight, on, nspike, tzero (ms),x,y, z, u, tsyn (ms)) {

	INITIAL {
		x = 0
		y = 0
		z = 0
		u = 0 :u0
		tsyn = t
		nspike = 1
	}

	if(onSET){on=0 onSET=0}
        if (flag == 0) { 
		: Qui faccio rientrare la modulazione presinaptica
		nspike = nspike + 1
		if (!on) {
			tzero = t	
			tpre=t	: activates diffusion
			on = 1				
			z = z*exp(-(t - tsyn)/tau_rec)		
			z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec)) / ((tau_1/tau_rec)-1) )
			y = y*exp(-(t - tsyn)/tau_1)			
			x = 1-y-z
				
			if (tau_facil > 0) { 
				u = u*exp(-(t - tsyn)/tau_facil)
				u = u + U * ( 1 - u )							
			} else { u = U }
			
			y = y + x * u
			T=Tmax*y
				
			PRE[numpulses]=y
			tspike[numpulses]=t
			numpulses=numpulses+1
			tsyn = t			
		}
		net_send(Cdur, nspike)
        }
	if (flag == nspike) { 
			T = 0
			on = 0
	}
}

