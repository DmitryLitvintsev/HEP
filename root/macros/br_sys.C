{
	double n1  = 36;
	double dn1 = 8;

	double n0  = 2846;
	double dn0 = 70;

	double n2  = 60.;
	double dn2 = 6.;

	double n3  = 24.7;
	double dn3 = 6.4;
	
	double epi_1  = 0.642;
	double depi_1 = 0.12;

	double epi_2 = 0.690;
	double depi_2 = 0.14;

	double epi_3 = 0.80;
	double depi_3 = 0.12;
	
	double r10 = 0.67;
	double dr10 = 0.06;
	
	double r20 = 0.64;
	double dr20 = 0.06;

	double r30 = 0.344;
	double dr30 = 0.03;

	double d = n0-n1/epi_1-n2/epi_2-3.*n3/epi_3;

	double r1 = n1/r10/epi_1/d;	

	double dr1_dn1  = 1./r10/epi_1/d+n1/r10/pow(epi_1,2)/d/d;
	double dr1_dr10 = - n1/epi_1/d/pow(r10,2);
	double dr1_de1  = - n1 / r10 / pow(epi_1,2) / d - n1 / r10 / pow(epi_1,3) * n1 / d / d;
	double dr1_de2  = - n1 / r10 / epi_1 / pow(epi_2,2) * n2 / d / d;
	double dr1_de3  = - n1 / r10 / epi_1 / pow(epi_3,2) * 3.*n3 / d / d;
	double dr1_dn0  = - n1 / r10 / epi_1 / d / d;
	double dr1_dn2 = n1 / r10 / epi_1 / epi_2 / d / d ;
	double dr1_dn3 = 3.* n1 / r10 / epi_1 / epi_3 / d / d ;

	double dr1 = sqrt(pow(dr1_dn1,2)*pow(dn1,2)+
			  pow(dr1_dr10,2)*pow(dr10,2)+
			  pow(dr1_de1,2)*pow(depi_1,2)+
			  pow(dr1_de2,2)*pow(depi_2,2)+
			  pow(dr1_de3,2)*pow(depi_3,2)+
			  pow(dr1_dn0,2)*pow(dn0,2)+
			  pow(dr1_dn2,2)*pow(dn2,2)+
			  pow(dr1_dn3,2)*pow(dn3,2));

// 2/3 - Br(lc* -> Lc 2pi)

	printf("B_1/B_0: %8.6e +- %8.6e\n",r1*3./2.,dr1*3./2.);
					  

	double r2 = n2/r20/epi_2/d;


	double dr2_dn2  = 1./r20/epi_2/d+n2/r20/pow(epi_2,2)/d/d;
	double dr2_dr20 = - n2/epi_2/d/pow(r20,2);
	double dr2_de2  = - n2 / r20 / pow(epi_2,2) / d - n2 / r20 / pow(epi_2,3) * n2 / d / d;
	double dr2_de1  = - n2 / r20 / epi_2 / pow(epi_1,2) * n1 / d / d;
	double dr2_de3  = - n2 / r20 / epi_2 / pow(epi_3,2) * 3 * n3 / d / d;
	double dr2_dn0  = - n2 / r20 / epi_2 / d / d;
	double dr2_dn1 = n2 / r20 / epi_2 / epi_1 / d / d ;
	double dr2_dn3 = 3. * n2 / r20 / epi_2 / epi_3 / d / d ;

	double dr2 = sqrt(pow(dr2_dn1,2)*pow(dn1,2)+
			  pow(dr2_dr20,2)*pow(dr20,2)+
			  pow(dr2_de1,2)*pow(depi_1,2)+
			  pow(dr2_de2,2)*pow(depi_2,2)+
			  pow(dr2_de3,2)*pow(depi_3,2)+
			  pow(dr2_dn0,2)*pow(dn0,2)+
			  pow(dr2_dn2,2)*pow(dn2,2)+
			  pow(dr2_dn3,2)*pow(dn3,2));


	printf("B_2/B_0: %8.6e +- %8.6e\n",r2*3./2.,dr2*3./2.);
					  

	double r3 = n3/r30/epi_3/d;

	double dr3_dn3  = 1./r30/epi_3/d+3.*n3/r30/pow(epi_3,2)/d/d;
	double dr3_dr30 = - n3/epi_3/d/pow(r30,2);
	double dr3_de3  = - n3 / r30 / pow(epi_3,2) / d - 3.* n3 / r30 / pow(epi_3,3) * n3 / d / d;
	double dr3_de1  = - n3 / r30 / epi_3 / pow(epi_1,2) * n1 / d / d;
	double dr3_de2  = - n3 / r30 / epi_3 / pow(epi_2,2) * n2 / d / d;
	double dr3_dn0  = - n3 / r30 / epi_3 / d / d;
	double dr3_dn1 = n3 / r30 / epi_3 / epi_1 / d / d ;
	double dr3_dn2 = n3 / r30 / epi_3 / epi_2 / d / d ;

	double dr3 = sqrt(pow(dr3_dn1,2)*pow(dn1,2)+
			  pow(dr3_dr30,2)*pow(dr30,2)+
			  pow(dr3_de1,2)*pow(depi_1,2)+
			  pow(dr3_de2,2)*pow(depi_2,2)+
			  pow(dr3_de3,2)*pow(depi_3,2)+
			  pow(dr3_dn0,2)*pow(dn0,2)+
			  pow(dr3_dn2,2)*pow(dn2,2)+
			  pow(dr3_dn3,2)*pow(dn3,2));

	printf("B_3/B_0: %8.6e +- %8.6e\n",r3,dr3);



        
}
