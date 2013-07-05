TITLE NMDA sinaptico

COMMENT
	NMDA model from article (version 15 sept 2004).
	Based on Nieus et al, 2006.
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDAS
	NONSPECIFIC_CURRENT i
	RANGE Rb,Ru,Rd,Rr,Ro, Rc,rb
	RANGE g,gmax,Cdur,Erev 			
	RANGE MgBlock,v0_block,k_block
	RANGE gg1,gg2,gg3
	RANGE T,Trelease,Tdiff
	RANGE tau_1,tau_rec,tau_facil,U,u0	
	RANGE A1,A2,A3,tau_dec1,tau_dec2,tau_dec3		: comes from fit 
	RANGE tdelay,ton,PRE
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
	PI	= (pi)		(1)
}

PARAMETER {
	Rb		=  5		(/ms/mM)  	: binding  
	Ru		=  0.1		(/ms)		: unbinding
	Rd		=  12e-4  	(/ms)		: desensitization
	Rr		=  9e-3		(/ms)		: resensitization 
	Ro		=  3e-2 	(/ms)		: opening
	Rc		=  0.966	(/ms)		: closing	
	Erev		= -3.7  	(mV)	: 0 (mV)
	gmax		= 10e4  	(pS)	: 7e3 : 4e4
	v0_block 	= -20 		(mV)	: -8.69 (mV)	: -18.69 (mV) : -32.7 (mV)
	k_block 	= 13		(mV)

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
	u0 		= 0 (1) 	< 0, 1 >	: se u0=0 al primo colpo y=U}

ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: actual conductance
	rb		(/ms)    : binding
	MgBlock

	T		(mM)
	Trelease	(mM)
	Tdiff		(mM)
	Tdiff_0		(mM)
	tdelay		(ms)
	ton		(ms)
	x
	PRE	
}

STATE {
	C0		: unbound
	C1		: single bound
	C2		: double bound
	D		: desensitized
	O		: open
	gg1
	gg2
	gg3
	sink
}

INITIAL {
	rates(v)
	C0		=	1
	C1		=	0
	C2		=	0
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
	rates(v)
	
	if( (t-ton)>tdelay) {
		Tdiff=gg1+gg2+gg3
		Tdiff_0 = Tdiff
	}else{
		Tdiff=Tdiff_0+(A1+A2+A3)*PRE*(t-ton)/tdelay
	}
	Trelease	= 	Tdiff 
	SOLVE kstates METHOD sparse
	g = gmax* O 			
	i = (1e-6) * g * (v - Erev) * MgBlock 
}

KINETIC kstates {	
	rb = Rb * Trelease 
	~ C0 <-> C1	(rb,Ru) 	
	~ C1 <-> C2	(rb,Ru)		
	~ C2 <-> D	(Rd,Rr)
	~ C2 <-> O	(Ro,Rc)
	CONSERVE C0+C1+C2+D+O = 1
	: Glutamate diffusion wave
	~ gg1 <-> sink (1/tau_dec1,0)
	~ gg2 <-> sink (1/tau_dec2,0)
	~ gg3 <-> sink (1/tau_dec3,0)
}

PROCEDURE rates(v(mV)) {
	: E' necessario includere DEPEND v0_block,k_block per aggiornare le tabelle!
	TABLE MgBlock DEPEND v0_block,k_block FROM -120 TO 30 WITH 150
	MgBlock = 1 / ( 1 + exp ( - ( v - v0_block ) / k_block ) )
}


NET_RECEIVE(weight, on, nspike,flagtdel, t0 (ms),y, z, u, tsyn (ms)) {
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
			tsyn=t
			PRE=y
		}
		net_send(tdelay,flagtdel)						
    	}
	if (flag == flagtdel){
		flagtdel = flagtdel+1
		state_discontinuity(gg1,gg1+A1*x*u)	 
		state_discontinuity(gg2,gg2+A2*x*u)	 
		state_discontinuity(gg3,gg3+A3*x*u)	 
		on=0
	}
}	 

 
