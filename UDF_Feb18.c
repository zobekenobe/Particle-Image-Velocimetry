#include "udf.h"
#include "sg.h"

real m_A0 = 0.0002; 	// solute mass fraction at inlet (kg/kg)
real Re = 400.; 			// feed Reynolds Number - dimensionaless
real h = 0.0025; 		// channel height 
real l = 0.060; 			// channel length 
real b = 0.030;	
real Dh = 4.0*(30.0/31.0); // Hydraulic Diameter
real J_v = -0.73e-05; 	// Permeate flux at the wall
real m_AP = 5.4e-05	;	// m_AP mass fraction of solute on permeate side

// user defined memory C_UDM
int x_nondim = 0; 		// x/h dimensionaless distance from the wall
int y_nondim = 1; 		// y/h dimensionaless distance fron the wall
int gamma    = 2;		// concentration polarization layer thickness (m)
int delta 	 = 3;		// concentration polarization layer thickness (eq 11)
int k_delta	 = 4; 		// mass transfer coefficient (eq 13)

// user defined scalar C_UDS index
int m_A      = 0;

// Initialize the domain
DEFINE_INIT(solute_init, domain_pointer)
{
	Thread *thread_pointer; 
	cell_t cellt;

	thread_loop_c(thread_pointer, domain_pointer)
	{
		begin_c_loop(cellt, thread_pointer)
		{
			C_UDSI(cellt, thread_pointer, m_A) = m_A0;
		}
		end_c_loop(cellt, thread_pointer)
	}
}

DEFINE_PROFILE(Side_Inlet_Velocity, thread_pointer, position_index)
{
	face_t facet;
	cell_t cellt;
	
	/*real rho = 997.1 + 694.*m_A0;
	real vis = 0.89e-03*(1+1.63*m_A0);*/

	real rho = 997.1 + 694.*C_UDSI(cellt, thread_pointer, m_A);
	real vis = 0.89e-03*(1+1.63*C_UDSI(cellt, thread_pointer, m_A));


	begin_f_loop(facet, thread_pointer)
	{
		F_PROFILE(facet, thread_pointer, position_index) = Re*vis/(rho*Dh);
	}
	end_f_loop(facet, thread_pointer)
}

DEFINE_PROFILE(Top_Inlet_Velocity, thread_pointer, position_index)
{
	face_t facet;
	cell_t cellt;
	
	/*real rho = 997.1 + 694.*m_A0;
	real vis = 0.89e-03*(1+1.63*m_A0);*/

	real rho = 997.1 + 694.*C_UDSI(cellt, thread_pointer, m_A);
	real vis = 0.89e-03*(1+1.63*C_UDSI(cellt, thread_pointer, m_A));

	begin_f_loop(facet, thread_pointer)
	{
		F_PROFILE(facet, thread_pointer, position_index) = (Re * vis )/(10 * rho*Dh);
	}
	end_f_loop(facet, thread_pointer)
}

DEFINE_PROFILE(Side_inlet_conc, thread_pointer, position_index)
{
	face_t facet;

	begin_f_loop(facet, thread_pointer)
	{
		F_PROFILE(facet, thread_pointer, position_index) = m_A0;
	}
	end_f_loop(facet, thread_pointer)
}

DEFINE_PROFILE(bottom_outlet, thread_pointer, position_index)
{
	// mixed boundary condition
	Thread *t0 = thread_pointer->t0;
	cell_t c0;
	face_t facet;

	real A[ND_ND];
	real ds; 			// distance between  cell and the face centroids
	real es[ND_ND];
	real A_by_es;
	real dr0[ND_ND];
	real source;

	real dphi[ND_ND];
	real k;

	real temp1, temp2; 
	real m_Aw; 		   // m_A at the wall

	begin_f_loop(facet, thread_pointer)
	{
		BOUNDARY_FACE_GEOMETRY(facet, thread_pointer, A, ds, es, A_by_es, dr0);

		c0 = F_C0(facet, thread_pointer);
		k = C_UDSI_DIFF(c0, t0, m_A);

		if(NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(m_A)))) source = 0.0;
			
		else
			BOUNDARY_SECONDARY_GRADIENT_SOURCE(source, SV_UDSI_G(m_A), dphi, es, A_by_es, k);

		temp1 = k*A_by_es/ds;
		temp2 = J_v * C_R(c0, t0) * NV_MAG(A);
		m_Aw  = (temp1 * C_UDSI(c0, t0, m_A) - source + temp2 * m_AP)/(temp1 + temp2);

		F_PROFILE(facet, thread_pointer, position_index) = m_Aw;

	}
	end_f_loop(facet,  thread_pointer)

}

/* Trying to find the mass flow rate at the bottom wall*/
DEFINE_UDS_FLUX(uds_flux, f, t, i)
{
	cell_t c0, c1 = -1;
	Thread *t0, *t1 = NULL;

	real NV_VEC(psi_vec), NV_VEC(A), flux = 0.0;

	c0 = F_C(f,t);
	t0 = F_C0_THREAD(f,t);
	F_AREA(A,f,t);

	/*if the face lies at the domain boundary, use the face values*/
	/*if the face lies in the domain, use the average of the cell values*/

	if(BOUNDARY_FACE_THREAD_P(t))
	{
		real dens;

		/*Depending on the boundary condition, density may not be defined on the boundary face thread*/
		if(NNULLP(THREAD_STORAGE(t, SV_DENSITY)))
			/*Density is set to the face value*/
			dens = F_R(f,t);
		else
			/*Density is set to the cell value*/
			dens = C_R(c0,t0);

		NV_DS(psi_vec, =,  F_U(f,t), F_V(f, t), F_W(f,t), *, dens);

		flux = NV_DOT(psi_vec, A); /*flux through the Face prescribed*/
	}
	else
	{
		c1 = F_C1(f,t);
		t1 = F_C1_THREAD(f, t);

		NV_DS(psi_vec, =, C_U(c0,t0), C_V(c0,t0),C_W(c0,t0),*,C_R(c0,t0));
		NV_DS(psi_vec, += C_U(c1,t1), C_V(c1,t1),C_W(c1,t1),*,C_R(c1,t1));

		flux = NV_DOT(psi_vec, A)/2.0; /*Average flux through the face*/

		return flux;
	}

}

/*Define the mass transfer at the interface surface*/
DEFINE_EXCHANGE_PROPERTY(interface_exchange, sourceIndex, sinkIndex, sourceSpeciesFrom, sinkS)
{

}

DEFINE_MASS_TRANSFER(interface_exchange, mixture_thread, fromPhaseIndex, fromSpeciesIndex, toPhaseindex, toSpeciesIndex)
{
		
}

/* defining the fluid density (page 2090) */

DEFINE_PROPERTY(density_solvent, cellt, thread_pointer)
{
	return 997.1 + 694.*C_UDSI(cellt, thread_pointer, m_A);
}

/* defining the solute diffusivity in the solvent */ 
DEFINE_DIFFUSIVITY(diffusivity, cellt, thread_pointer, position_index)
{
	return C_UDSI(cellt, thread_pointer, m_A) > 0.006 ? C_R(cellt, thread_pointer)*1.45e-09: C_R(cellt, thread_pointer)*1.61e-09*(1 - 14 * C_UDSI(cellt, thread_pointer, m_A));
}

DEFINE_PROPERTY(solvent_viscosity, cellt, thread_pointer)
{
	return 0.89e-03*(1+0.63*C_UDSI(cellt, thread_pointer, m_A));
}

DEFINE_ADJUST(performanceCalculation, domain_pointer)
{
	/*This routine calculates the x/h, y/h, gamma, delta and k/delta*/

	//Domain* domain_pointer;
	Thread* thread_pointer;
	cell_t cellt;

	real x[ND_ND]; // store the coordinates of the cell centroid

	//Domain* domain_pointer = Get_Domain(1);

	thread_loop_c(thread_pointer, domain_pointer)
	{
		begin_c_loop(cellt, thread_pointer)
		{
			C_CENTROID(x, cellt, thread_pointer);

			C_UDMI(cellt, thread_pointer, x_nondim) = x[0]/h;

			C_UDMI(cellt, thread_pointer, y_nondim) = x[1]/h;

			C_UDMI(cellt, thread_pointer, gamma) = C_UDSI(cellt, thread_pointer, m_A)/m_A0 - 1.0;

			C_UDMI(cellt, thread_pointer, delta) = log((C_UDSI(cellt, thread_pointer, m_A) - m_AP)/(m_A0 - m_AP))*C_UDSI_DIFF(cellt, thread_pointer, m_A)/(J_v * C_R(cellt, thread_pointer))/h;

			C_UDMI(cellt, thread_pointer, k_delta) = C_UDSI_DIFF(cellt, thread_pointer, m_A)/C_UDMI(cellt, thread_pointer, delta); 	
			
		}
		end_c_loop(cellt, thread_pointer)

	}
		//Message("This is working fine \n and this is the value of h %g", C_UDMI(cellt, thread_pointer, x_nondim));	

}

