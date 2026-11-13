TOPONIUM K-FACTORS vs TOP MASS, pp -> t tbar, 13 TeV
=====================================================

Contents
  mt_<m_t>/<distributionname>.csv   the 7 distributions of the 13 TeV run card, bin edges identical to the card:
      ATLAS_TTBAR_13TEV_HADR_DIF_YTTBAR  |y_ttbar|  0..2.4 (12 bins)
      CMS_TTBAR_13TEV_2L_DIF_YTTBAR      y_ttbar   -2.6..2.6 (10 bins)
      CMS_TTBAR_13TEV_LJ_DIF_MTTBAR      m_ttbar    250..3500 (15 bins)
      ATLAS_TTBAR_13TEV_LJ_DIF_MTTBAR    m_ttbar    325..2000 (9 bins)
      ATLAS_TTBAR_13TEV_LJ_DIF_YTTBAR    |y_ttbar|  0..2.5 (7 bins)
      ATLAS_TTBAR_13TEV_LJ_DIF_YBTTB     |y_boost| = |y_t + y_tbar|/2, 0..2.5 (9 bins)
      ATLAS_TTBAR_13TEV_LJ_DIF_HTTTB     H_T = pT(t) + pT(tbar), 0..2000 (9 bins)
  mt_<m_t>/other_binnings/          five further binnings (pT_top, |y_top|, 2D |y_top| x pT_top, m_ttbar 345-1600),
                                    not in the 13 TeV run card
  m_t = 170, 171.25, 172.5, 173.75, 175 GeV
  K_vs_mt_summary.txt               K (+- MC error) for every bin with an effect, one column per mass (old names:
                                    m_ttbar_250-3500 = CMS_..._MTTBAR, m_ttbar_325-2000 = ATLAS_..._MTTBAR, etc.)

How to use
  sigma_bin(with toponium) = K x sigma_bin(SM), bin by bin, at the m_t used for the SM prediction.
  K_norm = K / K(whole binning): use it for distributions normalised to their own integral.
  Above the threshold region K = 1 exactly (all m_ttbar bins >= 400 GeV).

CSV columns
  <obs>_lo, <obs>_hi   bin edges (lower included, upper excluded); 2D binning: |y_t| edges, then pT_t edges
  K                    1 + Delta/C (see below)
  stat_err             Monte Carlo error on K
  model_unc            |K(p* < 100 GeV) - K(p* < 50 GeV)|: sensitivity to how far in p* the non-relativistic
                       Green's function is trusted. A convention-based spread, not a 1-sigma uncertainty.
  K_norm               K / K over the whole binning

Set-up (identical at every mass apart from m_t, Gamma_t)
  Process     p p -> t tbar -> (b l+ nu)(b~ l- nu~), l = e, mu; LO MadGraph5_aMC@NLO 2.9.27; PDF CT25NNLO (its alpha_s).
              Tops reconstructed from their decay products (off-shell, parton level).
  Width       Gamma_t from the LO formula, normalised to 1.49 GeV at m_t = 173:
              170: 1.4011, 171.25: 1.4378, 172.5: 1.4750, 173.75: 1.5127, 175: 1.5509 GeV.
  Threshold   |M|^2 modified near threshold with Green's functions (Fuks et al., arXiv:2411.18962):
                gg:    |M_1|^2 |G1/G0|^2 + |M_8|^2 |G8/G0|^2   (colour singlet / octet parts)
                qqbar: |M|^2 |G8/G0|^2                         (pure octet at LO)
              in the region 2 m_t - 20 <= W <= 2 m_t + 4 GeV (W = m_ttbar) and p* < 80 GeV (p* = top momentum in the
              ttbar rest frame); |M|^2 unchanged elsewhere.
  Singlet     JKT potential (Jezabek, Kuhn, Teubner, Z. Phys. C 56 (1992) 653) with its own two-loop alpha_s,
              alpha_s(mZ) = 0.12, NLO bracket, JKT infrared form below 5 GeV.
  Octet       repulsive Coulomb, alpha_s fixed at CT25NNLO's alpha_s(25 GeV) = 0.1469.
  K           C = SM cross section in the bin from a full-phase-space LO sample (10M events per mass);
              Delta = (topo - SM) in gg + (topo - SM) in qqbar, from dedicated threshold-window samples
              (1M events each), counted in 2 m_t - 20 <= W <= 2 m_t + 4 and p* < 100 GeV.

Limitations (please read before using)
  1. LO throughout; K is a ratio of LO predictions.
  2. Extent of the threshold region (the "window").
     |M|^2 is modified only for 2 m_t - 20 <= W <= 2 m_t + 4 GeV (W = m_ttbar), i.e. from 20 GeV below to 4 GeV
     above threshold, and p* < 80 GeV: the convention of Fuks et al. (W <= 350 GeV at m_t = 173). Outside it,
     plain LO SM. The window moves with m_t, so every mass gets the same slice around its own threshold.
     For the run-card distributions:
       - first m_ttbar bins, CMS 250-400 and ATLAS 325-400: contain the whole window at every mass (320-344 GeV at
         m_t = 170 up to 330-354 GeV at 175), so their K is complete as defined and not affected by where the window
         edges fall. All higher m_ttbar bins have K = 1 exactly.
       - rapidity and H_T distributions: the window events spread over all bins; K varies smoothly, by < 0.001
         between 170 and 175 GeV.
     What the window leaves out: the Coulomb enhancement above E = +4 GeV (it falls off slowly; at m_t = 173 the
     last 2 GeV below the edge still carry ~15% of Delta). Extending the window is NOT obviously better: far above
     threshold the fixed-order SM prediction (NLO/NNLO) already contains the leading Coulomb terms, so modifying
     |M|^2 there risks double counting when K multiplies a higher-order SM prediction. The cut at +4 GeV keeps the
     correction to the region where resummation is needed (bound states, threshold), at the price of possibly
     missing some of the above-threshold effect. This is a matching choice, not a numerical error.
     (other_binnings/m_ttbar_345-1600: its 345 GeV edge lies inside the window, so its first bin's K mainly tracks
      how much of the window lies above 345 GeV, 1.0000 at m_t = 170 to 1.043 at 175. Not in the run card.)
  3. The singlet strength depends on the alpha_s normalisation: with alpha_s(mZ) = 0.12 (JKT) instead of 0.118,
     |G1/G0|^2 is 3-14% larger below threshold; K - 1 in the threshold bins is about 4% larger than with 0.118.
  4. model_unc covers only the p* variation; the items above are not included in it.

Numbers at a glance (K, MC error ~0.0001)
  bin                         m_t = 170   171.25   172.5   173.75   175
  CMS m_ttbar 250-400            1.0384   1.0404   1.0427   1.0452   1.0480
  ATLAS m_ttbar 325-400          1.0382   1.0404   1.0428   1.0452   1.0481
  ATLAS H_T 0-90                 1.0654   1.0656   1.0659   1.0661   1.0663
  ATLAS H_T 90-170               1.0014   1.0015   1.0016   1.0017   1.0018
  rapidity bins                  1.008 - 1.011 at every mass
  inclusive (all bins)  1.0088   1.0088   1.0087   1.0087   1.0086

Production (MSU HPCC, 2026-11-06): 70 MadGraph jobs (14 per mass), 70M events, 1 h 36 min wall time,
  295 core-hours; 9-34 min per 1M-event job on 16 cores.

Checks done on this run
  per mass: param_card m_t/Gamma_t, matrix-element threshold E = W - 2 m_t and window, generation window, G tables,
  14 samples (10x1M SM full phase space + 4x1M window), max W of window samples = 2 m_t + 4, qqbar beam-mirror
  symmetry (<y_ttbar> = 0 within errors), identical CSV layout (12 binnings, 110 bins) at every mass.
