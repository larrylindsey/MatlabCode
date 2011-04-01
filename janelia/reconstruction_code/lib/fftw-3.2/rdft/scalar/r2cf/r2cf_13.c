/*
 * Copyright (c) 2003, 2007-8 Matteo Frigo
 * Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Mon Nov 10 20:45:21 EST 2008 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_r2cf -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 13 -name r2cf_13 -include r2cf.h */

/*
 * This function contains 76 FP additions, 51 FP multiplications,
 * (or, 31 additions, 6 multiplications, 45 fused multiply/add),
 * 68 stack variables, 23 constants, and 26 memory accesses
 */
#include "r2cf.h"

static void r2cf_13(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP516520780, +0.516520780623489722840901288569017135705033622);
     DK(KP300462606, +0.300462606288665774426601772289207995520941381);
     DK(KP581704778, +0.581704778510515730456870384989698884939833902);
     DK(KP859542535, +0.859542535098774820163672132761689612766401925);
     DK(KP769338817, +0.769338817572980603471413688209101117038278899);
     DK(KP686558370, +0.686558370781754340655719594850823015421401653);
     DK(KP514918778, +0.514918778086315755491789696138117261566051239);
     DK(KP251768516, +0.251768516431883313623436926934233488546674281);
     DK(KP503537032, +0.503537032863766627246873853868466977093348562);
     DK(KP904176221, +0.904176221990848204433795481776887926501523162);
     DK(KP575140729, +0.575140729474003121368385547455453388461001608);
     DK(KP957805992, +0.957805992594665126462521754605754580515587217);
     DK(KP600477271, +0.600477271932665282925769253334763009352012849);
     DK(KP522026385, +0.522026385161275033714027226654165028300441940);
     DK(KP301479260, +0.301479260047709873958013540496673347309208464);
     DK(KP226109445, +0.226109445035782405468510155372505010481906348);
     DK(KP853480001, +0.853480001859823990758994934970528322872359049);
     DK(KP083333333, +0.083333333333333333333333333333333333333333333);
     DK(KP612264650, +0.612264650376756543746494474777125408779395514);
     DK(KP038632954, +0.038632954644348171955506895830342264440241080);
     DK(KP302775637, +0.302775637731994646559610633735247973125648287);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, R0 = R0 + ivs, R1 = R1 + ivs, Cr = Cr + ovs, Ci = Ci + ovs, MAKE_VOLATILE_STRIDE(rs), MAKE_VOLATILE_STRIDE(csr), MAKE_VOLATILE_STRIDE(csi)) {
	  E T15, T1a, T11, T17, T14, T1b;
	  {
	       E TN, TD, TV, TA, Tb, TZ, T12, TS, Tx, Tu, Ti, TU;
	       TN = R0[0];
	       {
		    E T3, TP, Th, TB, Tp, Te, Tm, TC, Tr, T6, T9, Ts;
		    {
			 E Tn, Tf, Tg, T1, T2;
			 T1 = R0[WS(rs, 4)];
			 T2 = R1[WS(rs, 2)];
			 Tn = R0[WS(rs, 6)];
			 Tf = R0[WS(rs, 5)];
			 Tg = R0[WS(rs, 2)];
			 T3 = T1 - T2;
			 TP = T1 + T2;
			 {
			      E Tk, To, Tc, Td;
			      Tk = R1[0];
			      Th = Tf - Tg;
			      To = Tf + Tg;
			      Tc = R1[WS(rs, 4)];
			      Td = R1[WS(rs, 1)];
			      {
				   E T4, Tl, T5, T7, T8;
				   T4 = R1[WS(rs, 5)];
				   TB = Tn + To;
				   Tp = FMS(KP500000000, To, Tn);
				   Tl = Td + Tc;
				   Te = Tc - Td;
				   T5 = R0[WS(rs, 3)];
				   T7 = R1[WS(rs, 3)];
				   T8 = R0[WS(rs, 1)];
				   Tm = FNMS(KP500000000, Tl, Tk);
				   TC = Tk + Tl;
				   Tr = T4 + T5;
				   T6 = T4 - T5;
				   T9 = T7 - T8;
				   Ts = T7 + T8;
			      }
			 }
		    }
		    {
			 E TO, Ta, Tt, TQ;
			 TD = TB - TC;
			 TO = TC + TB;
			 Ta = T6 + T9;
			 TV = T6 - T9;
			 Tt = Tr - Ts;
			 TQ = Tr + Ts;
			 {
			      E TX, Tq, TR, TY;
			      TX = Tm - Tp;
			      Tq = Tm + Tp;
			      TA = T3 + Ta;
			      Tb = FNMS(KP500000000, Ta, T3);
			      TR = TP + TQ;
			      TY = FNMS(KP500000000, TQ, TP);
			      TZ = TX + TY;
			      T12 = TX - TY;
			      T15 = TO - TR;
			      TS = TO + TR;
			      Tx = FNMS(KP866025403, Tt, Tq);
			      Tu = FMA(KP866025403, Tt, Tq);
			      Ti = Te + Th;
			      TU = Th - Te;
			 }
		    }
	       }
	       Cr[0] = TN + TS;
	       {
		    E Tw, Tj, T13, TW;
		    Tw = FNMS(KP866025403, Ti, Tb);
		    Tj = FMA(KP866025403, Ti, Tb);
		    T13 = TU - TV;
		    TW = TU + TV;
		    {
			 E TE, TI, Tv, TF, TG, Ty;
			 TE = FMA(KP302775637, TD, TA);
			 TI = FNMS(KP302775637, TA, TD);
			 Tv = FMA(KP038632954, Tu, Tj);
			 TF = FNMS(KP038632954, Tj, Tu);
			 TG = FNMS(KP612264650, Tw, Tx);
			 Ty = FMA(KP612264650, Tx, Tw);
			 {
			      E TT, Tz, TK, TH, TM, T10, TL, TJ;
			      TT = FNMS(KP083333333, TS, TN);
			      Tz = FNMS(KP853480001, Ty, Tv);
			      TK = FMA(KP853480001, Ty, Tv);
			      TH = FNMS(KP853480001, TG, TF);
			      TM = FMA(KP853480001, TG, TF);
			      T1a = FNMS(KP226109445, TW, TZ);
			      T10 = FMA(KP301479260, TZ, TW);
			      TL = FNMS(KP522026385, Tz, TE);
			      Ci[WS(csi, 1)] = KP600477271 * (FMA(KP957805992, TE, Tz));
			      TJ = FMA(KP522026385, TH, TI);
			      Ci[WS(csi, 5)] = -(KP600477271 * (FNMS(KP957805992, TI, TH)));
			      Ci[WS(csi, 4)] = -(KP575140729 * (FMA(KP904176221, TM, TL)));
			      Ci[WS(csi, 3)] = KP575140729 * (FNMS(KP904176221, TM, TL));
			      Ci[WS(csi, 6)] = KP575140729 * (FMA(KP904176221, TK, TJ));
			      Ci[WS(csi, 2)] = KP575140729 * (FNMS(KP904176221, TK, TJ));
			      T11 = FMA(KP503537032, T10, TT);
			      T17 = FNMS(KP251768516, T10, TT);
			 }
			 T14 = FNMS(KP514918778, T13, T12);
			 T1b = FMA(KP686558370, T12, T13);
		    }
	       }
	  }
	  {
	       E T1e, T1c, T18, T16, T1d, T19;
	       T1e = FMA(KP769338817, T1b, T1a);
	       T1c = FNMS(KP769338817, T1b, T1a);
	       T18 = FNMS(KP859542535, T14, T15);
	       T16 = FMA(KP581704778, T15, T14);
	       T1d = FNMS(KP300462606, T18, T17);
	       T19 = FMA(KP300462606, T18, T17);
	       Cr[WS(csr, 1)] = FMA(KP516520780, T16, T11);
	       Cr[WS(csr, 5)] = FNMS(KP516520780, T16, T11);
	       Cr[WS(csr, 2)] = FMA(KP503537032, T1e, T1d);
	       Cr[WS(csr, 6)] = FNMS(KP503537032, T1e, T1d);
	       Cr[WS(csr, 3)] = FMA(KP503537032, T1c, T19);
	       Cr[WS(csr, 4)] = FNMS(KP503537032, T1c, T19);
	  }
     }
}

static const kr2c_desc desc = { 13, "r2cf_13", {31, 6, 45, 0}, &GENUS };

void X(codelet_r2cf_13) (planner *p) {
     X(kr2c_register) (p, r2cf_13, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_r2cf -compact -variables 4 -pipeline-latency 4 -n 13 -name r2cf_13 -include r2cf.h */

/*
 * This function contains 76 FP additions, 34 FP multiplications,
 * (or, 57 additions, 15 multiplications, 19 fused multiply/add),
 * 55 stack variables, 20 constants, and 26 memory accesses
 */
#include "r2cf.h"

static void r2cf_13(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP083333333, +0.083333333333333333333333333333333333333333333);
     DK(KP075902986, +0.075902986037193865983102897245103540356428373);
     DK(KP251768516, +0.251768516431883313623436926934233488546674281);
     DK(KP503537032, +0.503537032863766627246873853868466977093348562);
     DK(KP113854479, +0.113854479055790798974654345867655310534642560);
     DK(KP265966249, +0.265966249214837287587521063842185948798330267);
     DK(KP387390585, +0.387390585467617292130675966426762851778775217);
     DK(KP300462606, +0.300462606288665774426601772289207995520941381);
     DK(KP132983124, +0.132983124607418643793760531921092974399165133);
     DK(KP258260390, +0.258260390311744861420450644284508567852516811);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     DK(KP1_732050807, +1.732050807568877293527446341505872366942805254);
     DK(KP300238635, +0.300238635966332641462884626667381504676006424);
     DK(KP011599105, +0.011599105605768290721655456654083252189827041);
     DK(KP156891391, +0.156891391051584611046832726756003269660212636);
     DK(KP256247671, +0.256247671582936600958684654061725059144125175);
     DK(KP174138601, +0.174138601152135905005660794929264742616964676);
     DK(KP575140729, +0.575140729474003121368385547455453388461001608);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, R0 = R0 + ivs, R1 = R1 + ivs, Cr = Cr + ovs, Ci = Ci + ovs, MAKE_VOLATILE_STRIDE(rs), MAKE_VOLATILE_STRIDE(csr), MAKE_VOLATILE_STRIDE(csi)) {
	  E T13, Tb, Tm, TW, TX, T14, TU, T10, Tz, TB, Tu, TC, TR, T11;
	  T13 = R0[0];
	  {
	       E Te, TO, Ta, Tv, To, T5, Tw, Tp, Th, Tr, Tk, Ts, Tl, TP, Tc;
	       E Td;
	       Tc = R0[WS(rs, 4)];
	       Td = R1[WS(rs, 2)];
	       Te = Tc - Td;
	       TO = Tc + Td;
	       {
		    E T6, T7, T8, T9;
		    T6 = R1[0];
		    T7 = R1[WS(rs, 1)];
		    T8 = R1[WS(rs, 4)];
		    T9 = T7 + T8;
		    Ta = T6 + T9;
		    Tv = T7 - T8;
		    To = FNMS(KP500000000, T9, T6);
	       }
	       {
		    E T1, T2, T3, T4;
		    T1 = R0[WS(rs, 6)];
		    T2 = R0[WS(rs, 5)];
		    T3 = R0[WS(rs, 2)];
		    T4 = T2 + T3;
		    T5 = T1 + T4;
		    Tw = T2 - T3;
		    Tp = FNMS(KP500000000, T4, T1);
	       }
	       {
		    E Tf, Tg, Ti, Tj;
		    Tf = R1[WS(rs, 5)];
		    Tg = R0[WS(rs, 3)];
		    Th = Tf - Tg;
		    Tr = Tf + Tg;
		    Ti = R1[WS(rs, 3)];
		    Tj = R0[WS(rs, 1)];
		    Tk = Ti - Tj;
		    Ts = Ti + Tj;
	       }
	       Tl = Th + Tk;
	       TP = Tr + Ts;
	       Tb = T5 - Ta;
	       Tm = Te + Tl;
	       TW = Ta + T5;
	       TX = TO + TP;
	       T14 = TW + TX;
	       {
		    E TS, TT, Tx, Ty;
		    TS = Tv + Tw;
		    TT = Th - Tk;
		    TU = TS - TT;
		    T10 = TS + TT;
		    Tx = KP866025403 * (Tv - Tw);
		    Ty = FNMS(KP500000000, Tl, Te);
		    Tz = Tx + Ty;
		    TB = Ty - Tx;
	       }
	       {
		    E Tq, Tt, TN, TQ;
		    Tq = To - Tp;
		    Tt = KP866025403 * (Tr - Ts);
		    Tu = Tq - Tt;
		    TC = Tq + Tt;
		    TN = To + Tp;
		    TQ = FNMS(KP500000000, TP, TO);
		    TR = TN - TQ;
		    T11 = TN + TQ;
	       }
	  }
	  Cr[0] = T13 + T14;
	  {
	       E Tn, TG, TE, TF, TJ, TM, TK, TL;
	       Tn = FNMS(KP174138601, Tm, KP575140729 * Tb);
	       TG = FMA(KP174138601, Tb, KP575140729 * Tm);
	       {
		    E TA, TD, TH, TI;
		    TA = FNMS(KP156891391, Tz, KP256247671 * Tu);
		    TD = FNMS(KP300238635, TC, KP011599105 * TB);
		    TE = TA + TD;
		    TF = KP1_732050807 * (TD - TA);
		    TH = FMA(KP300238635, TB, KP011599105 * TC);
		    TI = FMA(KP256247671, Tz, KP156891391 * Tu);
		    TJ = TH - TI;
		    TM = KP1_732050807 * (TI + TH);
	       }
	       Ci[WS(csi, 5)] = FMA(KP2_000000000, TE, Tn);
	       Ci[WS(csi, 1)] = FMA(KP2_000000000, TJ, TG);
	       TK = TG - TJ;
	       Ci[WS(csi, 4)] = TF - TK;
	       Ci[WS(csi, 3)] = TF + TK;
	       TL = Tn - TE;
	       Ci[WS(csi, 2)] = TL - TM;
	       Ci[WS(csi, 6)] = TL + TM;
	  }
	  {
	       E TZ, T1b, T19, T1e, T16, T1a, TV, TY, T1c, T1d;
	       TV = FNMS(KP132983124, TU, KP258260390 * TR);
	       TY = KP300462606 * (TW - TX);
	       TZ = FMA(KP2_000000000, TV, TY);
	       T1b = TY - TV;
	       {
		    E T17, T18, T12, T15;
		    T17 = FMA(KP387390585, TU, KP265966249 * TR);
		    T18 = FNMS(KP503537032, T11, KP113854479 * T10);
		    T19 = T17 - T18;
		    T1e = T17 + T18;
		    T12 = FMA(KP251768516, T10, KP075902986 * T11);
		    T15 = FNMS(KP083333333, T14, T13);
		    T16 = FMA(KP2_000000000, T12, T15);
		    T1a = T15 - T12;
	       }
	       Cr[WS(csr, 1)] = TZ + T16;
	       Cr[WS(csr, 5)] = T16 - TZ;
	       T1c = T1a - T1b;
	       Cr[WS(csr, 2)] = T19 + T1c;
	       Cr[WS(csr, 6)] = T1c - T19;
	       T1d = T1b + T1a;
	       Cr[WS(csr, 3)] = T1d - T1e;
	       Cr[WS(csr, 4)] = T1e + T1d;
	  }
     }
}

static const kr2c_desc desc = { 13, "r2cf_13", {57, 15, 19, 0}, &GENUS };

void X(codelet_r2cf_13) (planner *p) {
     X(kr2c_register) (p, r2cf_13, &desc);
}

#endif				/* HAVE_FMA */