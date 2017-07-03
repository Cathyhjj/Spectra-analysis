
--------------------------------------------------------------------------------
-- Define the crystal field term.
--------------------------------------------------------------------------------
Dq_3d     = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, {{4, 0, -14}, {4, 3, -2 * math.sqrt(70)}, {4, -3, 2 * math.sqrt(70)}})
Dsigma_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, {{2, 0, -7}})
Dtau_3d   = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, {{4, 0, -21}})

if H_cf == 1 then
    Dq_3d_i     = -0.050
    Dsigma_3d_i = -0.010
    Dtau_3d_i   = -0.090

    Dq_3d_f     = -0.045
    Dsigma_3d_f = -0.009
    Dtau_3d_f   = -0.081

    H_i = H_i
         + Dq_3d_i     * Dq_3d
         + Dsigma_3d_i * Dsigma_3d
         + Dtau_3d_i   * Dtau_3d

    H_f = H_f
         + Dq_3d_f     * Dq_3d
         + Dsigma_3d_f * Dsigma_3d
         + Dtau_3d_f   * Dtau_3d
end

--------------------------------------------------------------------------------
-- Define the 3d-4p hybridization term.
--------------------------------------------------------------------------------
if H_3d_4p_hybridization == 1 then
    N_4p = NewOperator('Number', NFermions, IndexUp_4p, IndexUp_4p, {1, 1, 1})
         + NewOperator('Number', NFermions, IndexDn_4p, IndexDn_4p, {1, 1, 1})

    if H_coulomb == 1 then
        F0_3d_4p = NewOperator('U', NFermions, IndexUp_4p, IndexDn_4p, IndexUp_3d, IndexDn_3d, {1, 0}, {0, 0})
        F2_3d_4p = NewOperator('U', NFermions, IndexUp_4p, IndexDn_4p, IndexUp_3d, IndexDn_3d, {0, 1}, {0, 0})
        G1_3d_4p = NewOperator('U', NFermions, IndexUp_4p, IndexDn_4p, IndexUp_3d, IndexDn_3d, {0, 0}, {1, 0})
        G3_3d_4p = NewOperator('U', NFermions, IndexUp_4p, IndexDn_4p, IndexUp_3d, IndexDn_3d, {0, 0}, {0, 1})

        Delta_3d_4p_i = 11.653 * 1.0
        e_3d_i= -(NElectrons_3d - 1) * U_3d_3d_i / 2
        e_4p_i=  (NElectrons_3d - 1) * U_3d_3d_i / 2 + Delta_3d_4p_i

        Delta_3d_4p_f = 12.924 * 1.0
        e_3d_f= -(NElectrons_3d - 1) * U_3d_3d_f / 2
        e_4p_f=  (NElectrons_3d - 1) * U_3d_3d_f / 2 + Delta_3d_4p_f

        F2_3d_4p_i = 2.289 * 0.6
        G1_3d_4p_i = 0.853 * 0.6
        G3_3d_4p_i = 0.739 * 0.6

        F2_3d_4p_f = 0.0 * 0.6
        G1_3d_4p_f = 0.0 * 0.6
        G3_3d_4p_f = 0.0 * 0.6

        H_i = H_i
            + F2_3d_4p_i * F2_3d_4p
            + G1_3d_4p_i * G1_3d_4p
            + G3_3d_4p_i * G3_3d_4p
            + e_3d_i * N_3d
            + e_4p_i * N_4p

        H_f = H_f
            + F2_3d_4p_f * F2_3d_4p
            + G1_3d_4p_f * G1_3d_4p
            + G3_3d_4p_f * G3_3d_4p
            + e_3d_f * N_3d
            + e_4p_f * N_4p
    end

    if H_soc == 1 then
        ldots_4p = NewOperator('ldots', NFermions, IndexUp_4p, IndexDn_4p)

        zeta_4p_i = 0.057 * 1.0

        zeta_4p_f = 0.0 * 1.0

        H_i = H_i
            + zeta_4p_i * ldots_4p

        H_f = H_f
            + zeta_4p_f * ldots_4p
    end


    Va1_t2g_3d_4p = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_4p, IndexDn_4p, {{1, 0, -math.sqrt(3 / 5)}, {3, 0, -7 / math.sqrt(15)}})
                  + NewOperator('CF', NFermions, IndexUp_4p, IndexDn_4p, IndexUp_3d, IndexDn_3d, {{1, 0, -math.sqrt(3 / 5)}, {3, 0, -7 / math.sqrt(15)}})

    Ve_eg_3d_4p = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_4p, IndexDn_4p, {{1, 0, math.sqrt(6 / 5)}, {3, 0, -14 / 3 * math.sqrt(2 / 15)}, {3, 3, -7 / 3 / math.sqrt(3)}, {3, -3, 7 / 3 / math.sqrt(3)}})
                + NewOperator('CF', NFermions, IndexUp_4p, IndexDn_4p, IndexUp_3d, IndexDn_3d, {{1, 0, math.sqrt(6 / 5)}, {3, 0, -14 / 3 * math.sqrt(2 / 15)}, {3, 3, -7 / 3 / math.sqrt(3)}, {3, -3, 7 / 3 / math.sqrt(3)}})

    Ve_t2g_3d_4p = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_4p, IndexDn_4p, {{1, 0, math.sqrt(3 / 5)}, {3, 0, -14 / 3 / math.sqrt(15)}, {3, 3, 7 / 3 * math.sqrt(2 / 3)}, {3, -3, -7 / 3 * math.sqrt(2 / 3)}})
                 + NewOperator('CF', NFermions, IndexUp_4p, IndexDn_4p, IndexUp_3d, IndexDn_3d, {{1, 0, math.sqrt(3 / 5)}, {3, 0, -14 / 3 / math.sqrt(15)}, {3, 3, 7 / 3 * math.sqrt(2 / 3)}, {3, -3, -7 / 3 * math.sqrt(2 / 3)}})

	Va1_t2g_3d_4p_i = 0.90 * 1.0
	Ve_eg_3d_4p_i   = 0.20 * 1.0
	Ve_t2g_3d_4p_i  = 1.65 * 1.0

	Va1_t2g_3d_4p_f = 0.90 * 1.0
	Ve_eg_3d_4p_f   = 0.20 * 1.0
	Ve_t2g_3d_4p_f  = 1.65 * 1.0

    H_i = H_i
        + Va1_t2g_3d_4p_i * Va1_t2g_3d_4p
        + Ve_eg_3d_4p_i   * Ve_eg_3d_4p
        + Ve_t2g_3d_4p_i  * Ve_t2g_3d_4p

    H_f = H_f
        + Va1_t2g_3d_4p_f * Va1_t2g_3d_4p
        + Ve_eg_3d_4p_f   * Ve_eg_3d_4p
        + Ve_t2g_3d_4p_f  * Ve_t2g_3d_4p
end
