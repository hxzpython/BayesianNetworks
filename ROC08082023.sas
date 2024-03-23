dm "log; clear;";
options ps=80 ls=90 pageno=1 nodate formdlim='-' notes threads;
title;
ODS HTML CLOSE;
ODS HTML;
/*************************************************************************************/
/* Date: 08-08-2023                                                                  */
/* Note: (1) Program the penalized conditional score equation for variable selection */
/*           with measurement error models. Using the method in Ma and Li (2009).    */
/*       (2) On 07-22-09 I revised CSlin(07-20-09).sas to use Newton-Raphson to solve*/
/*           the penalized estimating equation.                                      */
/*       (3) On 07-29-09 I revised CSlin(07-22-09).sas to clean things up.           */
/*       (4) On 07-30-09 I revised CSlin(07-29-09).sas to drop those nearly-zero beta*/
/*           during the Newton-Raphson iteration.                                    */
/*       (5) On 08-05-09 I add the cross-validation procedure to choose lambda.      */
/*       (6) On 08-05-09 I revised CSlin(07-30-09).sas to add GCV and BIC selectors. */
/*       (7) On 08-09-09 I revised CSlin(08-05-09).sas to add another BIC selector.  */
/*       (8) On 08-11-09 I revised CSlin(08-09-09).sas to use the "correct" Fisher   */
/*           information definition to calculate the df.                             */
/*       (9) On 10-27-09 I revised CSlin(08-11-09).sas to see if there is any trend  */
/*           in any thing when the reliability ratio, or mev, varies.                */
/*       (10) On 11-23-09 I revised linSIMEX(10-27-09).sas to assume different varian*/
/*            ce and different mev for different x's.                                */
/*       (11) On 11-28-09 I revised linSIMEX(11-23-09).sas to implement the idea of  */
/*            using SIMEX to get a better BIC.                                       */
/*       (12) On 02-25-10 I revised linSIMEXBIC(11-28-09).sas to add the oracle CSE. */
/*       (13) On 03-02-10 I revised linSIMEXBIC(02-25-10).sas to use more c points so*/
/*            that I can make SIMEX plot.                                            */
/*       (14) On 03-04-10 I revised this program to use nonlinear extrapolation.     */
/*       (15) On 04-15-10 I revised linSIMEXBIC(03-02-10).sas to use nlpfdd to       */
/*            approximate the hessian.                                               */
/*       (16) On 10-01-10 revised linSIMEXBIC(04-15-10).sas to check empirically if  */
/*            naive BIC depends on sigmau.                                           */
/*       (17) On 02-22-11 revised linBIC(10-01-10).sas to program the revised BIC.   */
/*       (18) On 03-08-11 revised CS_lin_RBIC(02-22-11).sas to compare 6 revised BICs*/
/*       (19) On 04-02-11 revised CS_lin_RBIC(03-08-11).sas to include a method in   */
/*            Liang and Li (2009).                                                   */
/*       (20) On 04-12-11 revised CS_lin_RBIC(04-02-11).sas to include four more BICs*/
/*            defined on note 04/11/11 no. 2.                                        */
/*       (21) On 04-23-11 revised CS_lin_RBIC(04-12-11).sas to add a revised BIC for */
/*            Liang and Li's method.                                                 */
/*       (22) On 08-01-11 revised CS_lin_RBIC(04-23-11).sas to delete some competing */
/*            methods no longer considered in VSCondscore(05-06-11).pdf.             */
/*       (23) On 01-25-12 revised CS_lin_RBIC_cor(08-01-11).sas to correct the cal-  */
/*            culation of df_lambda described on note 01/24/12 no. 1.                */
/*       (24) On 02-01-12 revised CS_lin_RBIC_cor(01-25-12).sas to just keep on RR at*/
/*            at time.                                                               */
/*       (25) On 02-07-12 revised CS_lin_RBIC_onerr(02-01-12).sas to avoid skipping  */
/*            a lambda value when one of the considered methods does not converge at */
/*            this lambda value. I should let the other methods try on this lambda too*/
/*       (26) On 02-27-12 revised CS_lin_RBIC_onerr(02-07-12).sas to keep CS part of */
/*            the code only.                                                         */
/*       (27) On 02-28-12 revised CS_onerr(02-27-12).sas to do LS-based methods only.*/
/*       (29) On 02-29-12 revised LS_onerr(02-28-12).sas to use another LSCV defined */
/*            on note 02/29/12 no. 6.                                                */
/*       (30) On 03-05-12 revised LSCV(02-29-12).sas to delete CV part of the code.  */
/*       (31) On 01-24-2015 revised LSonly(03-05-12).sas to work on networks problem.*/
/*       (32) On 02-11-2015 revised LS(01-24-15).sas to create the weights wij more  */
/*            carefully as described on note 02/11/2015 no. 7.                       */
/*       (33) On 02-23-2015 revised CR(02-11-2015).sas to monitor false/true discovery*/
/*       (34) On 03-30-2015 revised CR(02-23-2015).sas to use a different way to     */
/*            choose lambda.                                                         */
/*       (35) On 04-26-2015 revised CR(03-30-2015).sas to implement the IRWLS algori-*/
/*            thm outlined on 04-21-2015 notes.                                      */
/*       (40) On 05-10-2015 revised CR_IRWLS(04-26-2015).sas to delete the debugging */
/*            lines.                                                                 */
/*       (41) On 06-14-2015 revised CR_IRWLS(05-10-2015).sas to fix the "weighting"  */
/*            matrix in SIC from iteration to interation (motivated by Dziak (2006). */
/*       (42) On 06-30-2015 revised CR_IRWLS(06-14-2015).sas to use SCAD penalty.    */
/*       (43) On 07-06-2015 revised CR_SCAD(06-30-2015).sas to use the right updating*/
/*            formulas on 07/06/2015 notes.                                          */
/*       (44) On 07-24-2015 revised CR_SCAD(07-06-2015).sas to try SIC-type tuning   */
/*            parameter selector.                                                    */
/*       (45) On 07-31-2015 revised CR_SCADSIC(07-24-2015).sas to compare four dif-  */
/*            ferent SIC.                                                            */
/*       (46) On 08-19-2015 revised CR_SCADSIC4(07-31-2015).sas to try two more NPE- */
/*            based tuning parameter selectors.                                      */
/*       (47) On 08-31-2015 revised CR_SICNPE(08-19-2015).sas to drop two of the NPE-*/
/*            based tuning parameter selectors.                                      */
/*       (48) On 10-03-2015 revised CR_SCAD2(08-31-2015).sas to put back the SIC as  */
/*            that in Foyel & Drton (2010), as on 10-03-2015 no. 2 notes.            */ 
/*       (49) On 01-16-2016 revised CR_SCAD3(10-03-2015).sas to use unpenalized B in */
/*            model criteria.                                                        */
/*       (50) On 01-27-2016 revised CR_SCAD3(01-16-2016).sas to drop SIC4, the SIC in*/
/*            Foyel & Drton (2010), in order to focus on monitoring SIC1 and NPE.    */
/*       (51) On 02-06-2016 revised CR_SIC12(01-27-2016).sas to compare 4 tuning para*/
/*            meter selctors.                                                        */
/*       (52) On 03-21-2016 revised CR_SICNPE2(02-06-2016).sas to avoid repeatedly   */
/*            computing matrices that actually stay the same from iteration to iteration. */
/*       (53) On 04-13-2016 revised CR_SICNPE2(03-21-2016).sas to implementing testing*/
/*            differential networks.                                                 */
/*       (54) On 04-15-2016 revised diffnet_power(04-15-2016).sas to implement boots-*/
/*            trap.                                                                  */
/*       (55) On 04-23-2016 revised diffnet_powerboot(04-13-2016).sas to try another */
/*            test described on 04/13/2015 notes, pages 5~6.                         */
/*       (56) On 04-26-2016 revised diffnet(04-23-2016).sas to sum up the p indepen- */
/*            dent F statistics.                                                     */
/*       (57) On 04-29-2016 revised diffnet_sumF(04-26-2016).sas to use the same block*/
/*            of data to estimate a column of B and to evaluate a score vector.      */
/*       (58) On 05-01-2016 revised diffnet_sumF1(04-29-2016).sas to nj=k(p-1) such  */
/*            that the data structure used to compute the test statistic is parallel */
/*            to that used to estimate the graph structure.                          */
/*       (59) On 05-04-2016 revised diffnet_sumF(05-01-2016).sas to skip estimating  */
/*            the graph structure, and use the true graph structure under H0 in comput*/
/*            ing the test statistics.                                               */
/*       (60) On 05-05-2016 revised diffnet_sumFid(05-04-2016).sas to exclude the zero*/
/*            entries in the estimated B matrix when constructing the test statistics*/
/*       (61) On 06-09-2016 revised diffnet_sumFid(05-05-2016).sas to NOT avoid double*/
/*            dipping data in different scores, so that different scores are dependent*/
/*            even when the graph structure is not estimated. I did this because all */
/*            I did to avoid double dipping data in different scores, as done in     */
/*            diffnet_sumFid(05-05-2016).sas, still does not give me null distribution*/
/*            as the sum of p independent F statistics, and thus I have to use boots-*/
/*            trap to estimate the p-value anyway.                                   */
/*       (62) On 06-08-2016 revised diffnet_sumFid(06-08-2016).sas to perturb scores */
/*            in order to create bootstrap sample of the test statistic.             */
/*       (63) On 10-31-2016 revised diffnet_Gamidboot(06-08-2016).sas to add the code*/
/*            for graph structure estimation.                                        */
/*       (64) On 11-12-2016 revised diffnet_boot(10-31-2016).sas to implement another*/
/*            bootstrap idea that invovles resampling scores.                        */
/*       (65) On 01-18-2017 revised diffner_boot(11-12-2016).sas to use the data from*/
/*            graph2 to estimate B, given the structure estimated by data1. Also use */
/*            a different way to create bootstrap sample of the test stat based on   */
/*            data 1 with resampling applied to each block of interventional data.   */
/*       (66) On 01-25-2017 revised diffner_boot2(01-18-2017).sas to use the combined*/
/*            data for boostrapping.                                                 */
/*       (67) On 03-29-2017 revised diffnet_boot3(01-25-2017).sas to test structure  */
/*            same or different only.                                                */
/*       (68) On 05-06-2017 revised diffnet_struc(03-29-2017).sas to implement the   */
/*            wild bootstrap outlined on 05/06/2017.                                 */
/*       (69) On 06-02-2017 revised diffnet_struc(05-06-2017).sas to focus on no m.e.*/
/*            case first, so that generating bootstrap data is easier.               */
/*       (70) On 07-08-2017 revised diffnet_struc(06-02-2017).sas to transform and   */
/*            scale the esstimated residuals in the wild bootstrap.                  */
/*       (71) On 07-13-2017 revised diffnet_wild(07-08-2017).sas to try a different  */
/*            way to generate bootstrap data.                                        */
/*       (72) On 07-16-2017 revised diffnet_wild(07-13-2017).sas to combine wild boot*/
/*            trap and pair bootstrap in one program to save time in collecting simu-*/
/*            lation results.                                                        */
/*       (73) On 07-23-2017 revised wildpair(07-16-2017).sas to delete two extra     */
/*            steps of toposort before merging the graphs.                           */
/*       (74) On 08-13-2017 revised wildpair(07-23-2017).sas to try another wild boot*/
/*            strap method so that I can distinguish between two different sets of   */
/*            regression coefficients even though they have the same structure.      */
/*       (75) On 08-16-2017 revised wild2(08-13-2017).sas to try another wild boot-  */
/*            strap.                                                                 */
/*       (76) On 08-17-2017 revised wild2(08-16-2017).sas to use a slightlty different*/
/*            test statistic and wild bootstrap.                                     */
/*       (77) On 08-19-2017 revised wild2(08-17-2017).sas to use all data to estimate*/
/* 			  the graph to construct the test statistic.                             */
/*       (78) On 08-21-2017 revised wild2(08-17-2017).sas to generate bootstrap ver- */
/*            sion of X1 and X2.                                                     */
/*       (79) On 09-01-2017 revised wild2(08-21-2017).sas to change how to estimate  */
/*            model error variance used in the wild bootstrap.                       */
/*       (80) On 09-11-2017 revised wild2(09-11-2017).sas to use two lambda's to est-*/
/*            imate the graph structure using combined data, hoping to avoid underfit*/
/*            models.                                                                */
/*       (81) On 10-09-2017 revised wild2(09-11-2017).sas to implement yet another   */
/*            test stat described in notes 10-09-2017.                               */
/*       (82) On 05-28-2019 revised wild3(10-09-2017).sas to use the node-wise parent*/
/*            selection mathod and toplogical sorting based on pvalues to delete cyles*/
/*        	  implemented in Honest3_pv(03-26-2019).sas to estimate a network.       */
/*       (83) On 06-16-2019 revised diffnet(05-28-2019).sas to remove wild bootstrap */
/*            because I always get inflated Type I error from it. And I added Q1 with Q3. */
/*       (84) On 06-18-2019 revised diffnet13(06-16-2019).sas to put back Q2 and wild*/
/*            bootstrap.                                                             */
/*       (85) On 06-25-2019 revised diffnet123(06-18-2019).sas to revise Q1 and Q3 by*/
/*            not using the combined graph structure.                                */
/*       (86) On 06-27-2019 revised diffnet123(06-25-2019).sas to use a gamma distri-*/
/*            bution to approximate the null distribution of the test stat.          */
/*       (87) On 07-03-2019 revised diffnet123_gamma(06-27-2019).sas to try using a  */
/*            smaller nj in the wildbootstrap.                                       */
/*       (88) On 07-17-2019 revised diffnet123_nm(07-03-2019).sas to use two levels  */
/*            of smaller nj in the wildbootstrap in order to estimate the variance of*/
/*            the original test statistic.                                           */
/*       (89) On 07-18-2019 revised diffnet123_nm2(07-17-2019).sas to record wild    */
/*            bootstrap version of test statistics at three levels of nj.            */ 
/*       (90) On 07-24-2019 revised diffnet123_nm3(07-18-2019).sas to delete wild    */
/*            bootstrap vesion of test statistic based on nj-"2".                    */
/*       (91) On 02-19-2020 revised diffnet123_nmv(07-24-2019).sas to code the idea  */
/*            of identifying most influential nodes.                                 */
/*       (92) On 03-09-2020 revised diffner_dropnode(02-19-2020).sas to record how   */
/*            often other nodes rank top 1 and top 2, besides recording such frequen-*/
/*            cy for the two nodes I manipulate.                                     */
/*       (93) On 03-11-2020 revised findhub03092020.sas to have 100 MC rep per graph.*/
/*       (94) On 03-26-2020 revised findhub03112020.sas to NOT re-estimate the graph */
/*            structure when I delete a node.                                        */
/*       (99) On 05-19-2020 revised findhub03262020.sas to try the delete-one-node   */
/*            at a time idea but always use the true graph structure.                */
/*       (100) On 06-16-2020 revised findimpdel05192020.sas to try out the idea based*/
/*             on causal inference (Peters, Buhlman, and Meinshausen, 2016).         */
/*       (101) On 12-16-2020 revised findimpdel06162020.sas to use the average score */
/*             across several random splits to create two sets of experimental condi-*/
/*             tions.                                                                */
/*       (102) On 03-172021 revised findimp12162020.sas to combine this code with the*/
/*             other named findhub03262020.sas to implete all three methods we tried.*/
/*       (103) On 04-08-2021 revised DiNA03172021.sas to simulate two networks with  */
/*             user-specified driver nodes.                                          */
/*       (104) On 04-28-2021 revised DiNA04082021.sas to drop scores based on Q1 and */
/*             Q3.                                                                   */
/*       (105) On 04-28-2021 revised DiNA04282021.sas to allow two types of driver   */
/*             nodes.                                                                */
/*       (106) On 05-03-2021 revised DiNA04292021.sas to make a cleaning code with 3 */
/*             metrics to asssess the performance of each of 5 scores.               */
/*       (107) On 07-19-2021 revised DiNA05032021.sas to try a different PI score    */
/*             that combines the two PI scores by switching G1 and G2.               */
/*       (108) On 07-27-2021 revised DiNA07102021.sas to incorporate a normalization */
/*             factor in the PI scores.                                              */
/*       (109) On 09-13-2021 revised DiNA07272021.sas to generate G1 and G2 systema- */
/*             tically to control differential nodes and driver nodes.               */
/*       (110) On 10-11-2021 revised DiNA09132021.sas to use pvalue to decide on pick*/
/*             -ing top how many nodes based on PI scores.                           */
/*       (111) On 02-23-2022 revised DinNA10112021.sas to add DICERN method for comp-*/
/*             arison.                                                               */
/*       (112) On 03-04-2033 revised DiNADISCERN02232022.sas to add more threshold   */
/*             levels to decide the number of differential nodes, and also add one   */
/*             more competing method besides DISCERN.                                */
/*       (113) On 03-30-2022 revised DiNADISCERN03042022.sas to refine my mean test (*/
/*             based on the variance test result, and also to implement the permuta- */
/*  		   tion test for DISCERN.                                                */
/*       (114) On 04-01-2022 revised DiNA03302022.sas add permutation test for the   */
/*             competing DISCERN method.                                             */
/*       (115) On 04-04-2022 revised Dina04012022.sas to apply BH procedure on my p- */
/*  		   values from my PI method.                                             */
/*       (116) on 04-13-2022 revised DiNA04132022.sas to clean up the simulation and */
/*             record driver nodes identification from competing methods too.        */
/*       (117) On 06-21-2022 revised DiNA04132022.sas to see if re-estimating graph  */
/*             structures based on permutation data are needed.                      */
/*       (118) On 09-12-2022 revised DiNAperm06212022.sas to try parametric bootstrap*/
/*             to assess significance of PI and DISC.                                */
/*       (119) On 11-07-2022 revised DiNApermboot09122022.sas to use the combined data*/
/*             to get a common graph structure estimate.                             */
/*       (120) On 11-08-2022 revised DiNApermboot11072022.sas to use entries from the*/
/*             original Bhat matrix to fill in B matrix used to generate bootstrap data*/
/*       (121) On 11-30-2022 revised DiNApermboot11082022.sas to use resampling type */
/*             of bootstrap to get a power-like probability for each node.           */
/*       (122) On 08-08-2023 revised DiNApower11302022.sas to get ROC and AUC compar-*/
/*             ing PI and DISC.                                                      */
/*************************************************************************************/
/*libname out "G:\Research_Networks\Output";*/
libname out "C:\Research_Networks\Output";
/*libname out "E:\Research_Networks";*/
/*libname out "C:\Documents and Settings\huang\Research_Networks";*/

/*libname out '/work/huang24/Research_Network/sasoutput';*/


/*proc printto log='C:\Research_Networks\Output\mylog.txt';*/
/*run;*/

/*%put Site: &syssite Release: &sysvlong System: &sysscp &sysscpl;*/



data track;
    t=time(); d=today();
    put 'Time: ' t time9.;
    put 'Date: ' d date9.;
run;

proc iml;
start lsobj(beta, y, x) global(small3);
	n=nrow(y); 
   	score=x`*(y-x*beta)/n; 
	hh=0; 
	do i=1 to n;
		lof=x[i,]`*(y[i]-x[i,]*beta)-score; 
		hh=hh+lof*lof`;
    end;
	if (min(vecdiag(hh))>1e8) then do;
		sic=1e8; goto zou;
	end;
	sic=score`*sweep(hh/n)*score; 
	zou:
	return(sic);
finish; 

/* Estimate bigB for a sequence of candidate lambda's */ 
start Tuning1(lambda1, bestB1, theorder, lambmax, m, q, allx1, bigBini1) global(NoNode, nj);
	oldnj=nj; nj=nrow(allx1)/NoNode; 
	NoObsInt=nrow(allx1); NoObs=NoObsInt-nj; 
	bicpen1=log(NoObs)/NoObs;
	lambda1=lambmax;
	do mm=1 to m;
		lambda=lambmax*q**(mm-1); 
		
		bigBnow1=findgraph2(lambda, allx1, bigBini1); 
		pvmat=getpv2(bigBnow1, allx1);
		wantall=topsort(bigBnow1, pvmat);
		bigBnow1=wantall[, 1:NoNode];
		theorder=wantall[, (NoNode+1)];
	
		/* using the estimated graph, compute SIC tuning parameter selector (eveluated at the penalized Bhat). */
		sic1=0; 
		do j=1 to NoNode; 
			ignorej=loc((1:NoNode)^=j); 
			if j=1 then whichrow=(nj+1):NoObsInt;
			if (j>1 & j<NoNode) then whichrow=(1:((j-1)#nj))||((j#nj+1):NoObsInt);
			if j=NoNode then whichrow=1:((NoNode-1)#nj);
			
			yw=allx1[whichrow, ]; 
			yw=yw-yw[:,];
			y=yw[, j]; w=yw[, ignorej];
			sic1=sic1+lsobj(bigBnow1[ignorej, j], y, w);
		end;

		if mm=1 then do; 
			smallestsic1=sic1+sum(bigBnow1^=0)#bicpen1; 
			sic1=smallestsic1;
			bestB1=bigBnow1;
		end; 
		else do; 
			sic1=sic1+sum(bigBnow1^=0)#bicpen1; 
			if sic1<=smallestsic1 then do;
				smallestsic1=sic1;
				lambda1=lambda;
				bestB1=bigBnow1;
			end;
		end;
    end;

	nj=oldnj; 
finish;


/* Estimate bigB for a sequence of candidate lambda's */  
start Tuning(lambda1, lambda2, bestB1, bestB2, lambmax, m, q, allx1, allx2, bigBini1, bigBini2) global(NoNode, nj, NoObsInt, bicpen1, bicpen2);
	lambda1=lambmax; lambda2=lambmax;
	do mm=1 to m;
		lambda=lambmax*q**(mm-1); 
		
		bigBnow1=findgraph2(lambda, allx1, bigBini1); 
		pvmat=getpv2(bigBnow1, allx1);
		bigBnow1=topsort(bigBnow1, pvmat)[, 1:NoNode];
		
		bigBnow2=findgraph2(lambda, allx2, bigBini2); 
		pvmat=getpv2(bigBnow2, allx2);
		bigBnow2=topsort(bigBnow2, pvmat)[, 1:NoNode];
		
		/* using the estimated graph, compute SIC tuning parameter selector (eveluated at the penalized Bhat). */
		sic1=0; sic2=0; 
		do j=1 to NoNode; 
			ignorej=loc((1:NoNode)^=j); 
			if j=1 then whichrow=(nj+1):NoObsInt;
			if (j>1 & j<NoNode) then whichrow=(1:((j-1)#nj))||((j#nj+1):NoObsInt);
			if j=NoNode then whichrow=1:((NoNode-1)#nj);
			
			yw=allx1[whichrow, ]; 
			yw=yw-yw[:,];
			y=yw[, j]; w=yw[, ignorej];
			sic1=sic1+lsobj(bigBnow1[ignorej, j], y, w);

			yw=allx2[whichrow, ]; 
			yw=yw-yw[:,];
			y=yw[, j]; w=yw[, ignorej];
			sic2=sic2+lsobj(bigBnow2[ignorej, j], y, w);
		end;

		if mm=1 then do; 
			smallestsic1=sic1+sum(bigBnow1^=0)#bicpen1; 
			sic1=smallestsic1;
			bestB1=bigBnow1;

			smallestsic2=sic2+sum(bigBnow2^=0)#bicpen2; 
			sic2=smallestsic2;
			bestB2=bigBnow2;
		end; 
		else do; 
			sic1=sic1+sum(bigBnow1^=0)#bicpen1; 
			if sic1<=smallestsic1 then do;
				smallestsic1=sic1;
				lambda1=lambda;
				bestB1=bigBnow1;
			end;

			sic2=sic2+sum(bigBnow2^=0)#bicpen2; 
			if sic2<=smallestsic2 then do;
				smallestsic2=sic2;
				lambda2=lambda;
				bestB2=bigBnow2;
			end;
		end;
    end;
finish;


/* Use my earlier variable selection code to select parents for each node based on corrected score method with LQA algorithm*/
start findgraph2(lambda, allw, bigBini0) global(nj, threshold, maxit, small, small2, scada);
	pall=ncol(allw); dall=pall-1; totn=nrow(allw); 
	oldnj=nj; nj=totn/pall;
	nall=totn-nj; 
	 
	bigBnow=bigBini0; n=nall; itno=0; absdiff=10; 
	
	do while (absdiff>threshold & itno<=maxit);
		bigBlast=bigBnow; itno=itno+1;
 
		do j=1 to pall;
			ignorej=loc((1:pall)^=j);
			betanow=bigBnow[, j];
			
			if max(abs(betanow))=0 then do;
				rc=-1; goto doneold; 
			end; 

			if j=1 then whichrow=(nj+1):totn; 
			if (j>1 & j<pall) then whichrow=(1:((j-1)#nj))||((j#nj+1):totn);
			if j=pall then whichrow=1:(dall#nj);
			y=allw[whichrow, j]; y=y-y[:];
			
			label=loc(betanow^=0);
			betanow=betanow[label];
			d=ncol(label);
			
			w=allw[whichrow, label]; w=w-w[:,];

			absdiff_in=10; itno_in=0; rc=0;
		    do while (absdiff_in>threshold & itno_in<maxit);
				noneg=j(d, 1, 0);
				diff=scada#lambda-abs(betanow);
			    pick=loc(diff>0);
			    if ncol(pick)>0 then noneg[pick]=diff[pick];
				
			    penal=lambda#sign(betanow)#((abs(betanow)<=lambda)+(abs(betanow)>lambda)#noneg/((scada-1)#lambda));
			    siglam=diag(penal/abs(betanow));
				hess=-w`*w/n-siglam;
			    score=w`*(y-w*betanow)/n-penal;
			    
		        if (abs(det(hess))>threshold) then do;
		            parmest=betanow-inv(hess)*score;
		        	absdiff_in=max(abs(parmest-betanow)); itno_in=itno_in+1;
		            
		            drop=loc(abs(parmest)<small);
		            notdrop=loc(abs(parmest)>=small);

		            if (ncol(drop)>0 & ncol(notdrop)>0) then do;
		                parmest=parmest[notdrop];
		                w=w[, notdrop];
						label=label[notdrop];
		                d=nrow(label);
					end;

		            if ncol(notdrop)=0 then do;
		                rc=-1;  /* no predictor chosen */
		                goto doneold;
		            end;

					betanow=parmest;
		        end;
		        else do;
		            rc=-2; 
					goto doneold; /* The estimates in adjacent two steps are not close enough but Hessian is nearly singular. */
		        end;
			end;
		    doneold: 
			if rc=-1 then bigBnow[, j]=0;
			else do;
				bigBnow[label, j]=betanow;
				bigBnow[setdif(1:pall, label), j]=0;
			end;
		end;
		absdiff=max(abs(bigBlast-bigBnow));
	end;
	nj=oldnj;
	return(bigBnow); 
finish;

/* obtain a p-value matrix for the unpenalized estimated regression coefficients obtained from corrected score method */
start getpv2(graphest, allw) global(nj);
	pall=ncol(allw); totn=nrow(allw); 
	oldnj=nj; nj=totn/pall; 
	nall=totn-nj;

	pvmat=j(pall, pall, -8); 
	do j=1 to pall;
		whichcol=loc(graphest[, j]^=0); d=ncol(whichcol); 
		if (d>0) then do;
			if (j>1 & j<pall) then whichrow=(1:((j-1)#nj))||((j#nj+1):totn); 
			if j=1 then whichrow=(nj+1):totn;
			if j=pall then whichrow=1:((pall-1)#nj);
		
			/*~~~~~~~~~~~~ center but not standardize w and y ~~~~~~~~~~~~*/
			yw=allw[whichrow, j]||allw[whichrow, whichcol];
			yw=yw-yw[:,]; 
			y=yw[, 1]; w=yw[, 2:(d+1)]; 

			Bhat=inv(w`*w)*w`*y; 

			meatmat=0; breadmat=0; 
			do i=1 to nall;  
			    score=w[i, ]`*(y[i]-w[i,]*Bhat);
				hess=-w[i,]`*w[i,];
			    meatmat=meatmat+score*score`;
				breadmat=breadmat+hess;
			end;
			meatmat=meatmat/nall; breadmat=breadmat/nall;
			ainv=ginv(breadmat);

			varbhat=ainv*meatmat*ainv`/nall; 
			tval=Bhat/sqrt(vecdiag(varbhat));
			pvmat[whichcol, j]=2*CDF('T', -abs(tval), nall-d);
		end;		
	end;
	nj=oldnj;
	return(pvmat);  
finish;


/* delete cycles based on p-values only */
start topsort(graphest, pvmat); 
	if (norm(graphest)=0) then goto emptygraph;

/*	fudge=range(pvmat[loc(graphest^=0)])/1000;*/
	/* Find the ordering of nodes indicated by graphest, and if this ordering finding fails, it means graphest is not DAG */
	/* in which case, delete the entry of graphest that has the largest p-value and find the ordering of the resultant graph again. */
	
	NoNode=ncol(graphest);
	noid=(1:NoNode); topo=0;


	do while (ncol(topo)=1);
		anyrt=0;
		do j=1 to NoNode;
			if (norm(graphest[, j], "LInf")=0) then do;
				anyrt=anyrt+1;
				topo=topo||j; noid=noid[loc(noid^=j)]; /* one problem here: loc(noid^=j) may not exist, noid is 1 by 1 */
			end;	
		end;
		
		if (anyrt=0) then do;
			/* set the "weakest" link at 0 */
			weakest=loc(pvmat=max(pvmat[loc(graphest^=0)]));
			graphest[weakest]=0;
			pvmat[weakest]=-1;
		end;
	end;
	
	topo=topo[2:ncol(topo)]`; remnd=nrow(noid); 
	if (remnd=1) then topo=topo||noid;
	temp=graphest;
	do while (remnd>1);
		temp=graphest[noid, noid]; 
		temp_noid=noid;
		anyrt=0;
		do j=1 to remnd;
			thisnode=temp_noid[j];
			if (norm(temp[, j], "LInf")=0) then do;
				topo=topo||thisnode; 
				anyrt=anyrt+1;
				still=loc(noid^=thisnode);
				if (ncol(still)>0) then noid=noid[still];
			end;	
		end;
		if (anyrt=0) then do;
			/* set the "weakest" link at 0 */
			weakest=loc(pvmat=max(pvmat[loc(graphest^=0)]));
			graphest[weakest]=0;
			pvmat[weakest]=-1;
		end;
		remnd=nrow(noid);
		if (remnd=1 & ncol(topo)<NoNode) then topo=topo||noid;
	end;

	emptygraph:
	return(graphest||topo`);
finish;

/* delete cycles based on entries of regression coefficients in an initial regression coefficients matrix, return a DAG and a topological order */
start Khan(graphini); 
	NoNode=ncol(graphini);
	noid=(1:NoNode); topo=0;

	pvmat=abs(graphini); 

	do while (ncol(topo)=1);
		anyrt=0;
		do j=1 to NoNode;
			if (norm(graphini[, j], "LInf")=0) then do;
				anyrt=anyrt+1;
				topo=topo||j; noid=noid[loc(noid^=j)]; /* one problem here: loc(noid^=j) may not exist, noid is 1 by 1 */
			end;	
		end;
		
		if (anyrt=0) then do;
			/* set the "weakest" link at 0 */
			weakest=loc(pvmat=min(pvmat[loc(graphini^=0)]));
			graphini[weakest]=0;
			pvmat[weakest]=888;
		end;
	end;
	
	topo=topo[2:ncol(topo)]`; remnd=nrow(noid); 
	if (remnd=1) then topo=topo||noid;
	temp=graphini;
	do while (remnd>1);
		temp=graphini[noid, noid]; 
		temp_noid=noid;
		anyrt=0;
		do j=1 to remnd;
			thisnode=temp_noid[j];
			if (norm(temp[, j], "LInf")=0) then do;
				topo=topo||thisnode; 
				anyrt=anyrt+1;
				still=loc(noid^=thisnode);
				if (ncol(still)>0) then noid=noid[still];
			end;	
		end;
		if (anyrt=0) then do;
			/* set the "weakest" link at 0 */
			weakest=loc(pvmat=min(pvmat[loc(graphini^=0)]));
			graphini[weakest]=0;
			pvmat[weakest]=888;
		end;
		remnd=nrow(noid);
		if (remnd=1 & ncol(topo)<NoNode) then topo=topo||noid;
	end;

	
	return(graphini||topo`);
finish;

/* compute some metrics to assess the quality of an estimated graph */
start quality(graphest, truegraph);
	wanttotal=ncol(loc(truegraph^=0)); 
	pall=ncol(truegraph); 	
	allpossible=pall*(pall-1)/2; 
	allneg=allpossible-wanttotal;

	predicted=sum(graphest^=0);	/* number of edges in the estimated graph */
	if predicted=0 then do; /* this means the estimated graph is empty */	
		expected=0; reversed=0;
	end;
	else do; 
		expected=sum(truegraph[loc(graphest^=0)]^=0); 
		reversed=sum(truegraph`[loc(graphest^=0)]^=0);
	end;
	expected=sum(truegraph[loc(graphest^=0)]^=0); /* among the edges in the estimated graph, how many are also in the true graph */
	reversed=sum(truegraph`[loc(graphest^=0)]^=0); /* among the edges in the estimated graph, how many are reversed in the true graph */
	fp=predicted-expected-reversed; /* among the edges in the estimated graph, how many are absence in the true graph regardless of the direction */
	missing=sum((graphest[loc(truegraph^=0)]=0)#(graphest`[loc(truegraph^=0)]=0)); /* number of edges in the true graph that are missing in the estimated graph regardless of the direction */
	tpr=expected/wanttotal; /* true positive rate */
	if predicted=0 then fdr=0;
	else fdr=(reversed+fp)/predicted; /* false discovery rate */
	specif=(sum(graphest[loc(truegraph=0)]=0)-pall)/(sum(truegraph=0)-pall); /* among the zero off diag spots in the true B matrix, the proportion of spots that are also zero in the estimated B */
	corre=(expected+allneg*specif)/allpossible;  /* correctness rate */

	assess=tpr||fdr||specif||corre;
	return(assess);
finish;

/* compute pvalues used in PI score given two data set */
start meanvartest(sample1, sample2, samsize1, samsize2) global(varpvct); 
	sampvar1=var(sample1); sampvar2=var(sample2);
	
	larger=sampvar1<>sampvar2; 
	smaller=sampvar1><sampvar2;
	if larger>sampvar2 then do; 
		nlarge=samsize1; nsmall=samsize2;
	end; 
	else do; 
		nlarge=samsize2; nsmall=samsize1;
	end; 
	pvar=cdf("F", smaller/larger, nsmall-1, nlarge-1)+sdf("F", larger/smaller, nlarge-1, nsmall-1);

	if (pvar<varpvct) then do; /* two-sample t test without assuming equal variance */
		sumvar=sampvar1/samsize1+sampvar2/samsize2;
		df=sumvar**2/((sampvar1/samsize1)**2/(samsize1-1)+(sampvar2/samsize2)**2/(samsize2-1));
		tstat=(sample1[:]-sample2[:])/sqrt(sumvar);
	end; 
	else do; /* two-sample t test assuming equal variance */
		df=samsize1+samsize2-2; 
		poolvar=((samsize1-1)*sampvar1+(samsize2-1)*sampvar2)/df;
		tstat=(sample1[:]-sample2[:])/sqrt(poolvar*(1/samsize1+1/samsize2));
	end;	
	pmean=2*sdf("T", abs(tstat), df);

	lower=pmean><pvar; 

	return(pmean||pvar||lower);
finish; 

/* compute pvalues used in PI score given two data set but I do not do the mean test is the vairance test is significant */
start meanvartest2(sample1, sample2, samsize1, samsize2) global(varpvct); 
	sampvar1=var(sample1); sampvar2=var(sample2);
	
	larger=sampvar1<>sampvar2; 
	smaller=sampvar1><sampvar2;
	if larger>sampvar2 then do; 
		nlarge=samsize1; nsmall=samsize2;
	end; 
	else do; 
		nlarge=samsize2; nsmall=samsize1;
	end; 
	pvar=cdf("F", smaller/larger, nsmall-1, nlarge-1)+sdf("F", larger/smaller, nlarge-1, nsmall-1);

	if (pvar<varpvct) then pmean=-8; 
	else do; /* two-sample t test assuming equal variance */
		df=samsize1+samsize2-2; 
		poolvar=((samsize1-1)*sampvar1+(samsize2-1)*sampvar2)/df;
		tstat=(sample1[:]-sample2[:])/sqrt(poolvar*(1/samsize1+1/samsize2));
	end;	
	pmean=2*sdf("T", abs(tstat), df);

	return(pmean||pvar);
finish; 


/* Find PI scores and DISC scores of p nodes given two data sets from two networks and two estimated regression coefficients matrices */
start Find2score(PIscore, disc, tophowmany, allx1, allx2, bigBnow1, bigBnow2, randomseed) global(dall, NoNode, split, nj, NoObs, NoObsInt, pvcutoff); 
		/* Use the estimated graph structure to get unpenalized estimated regression coefficients */
		resids1to2=j(NoObs, NoNode, 0); resids2to1=j(NoObs, NoNode, 0); 
		resids1to1=j(NoObs, NoNode, 0); resids2to2=j(NoObs, NoNode, 0); 
		err12=j(NoObs, NoNode, 0); err21=j(NoObs, NoNode, 0); 
		do j=1 to NoNode;
			if (j>1 & j<NoNode) then whichrow=(1:((j-1)#nj))||((j#nj+1):NoObsInt); 
			if j=1 then whichrow=(nj+1):NoObsInt;
			if j=NoNode then whichrow=1:((NoNode-1)#nj);
			
			/* first assuming G1 structure */
			whichcol=loc(bigBnow1[, j]^=0); 

			d=ncol(whichcol);
			if (d>0) then do;
				/*~~~~~~~~~~~~ center but not standardize w and y ~~~~~~~~~~~~*/
				yw=allx2[whichrow, j]||allx2[whichrow, whichcol];
				yw=yw-yw[:,]; 
				y=yw[, 1]; w=yw[, 2:(d+1)];
				trouble1=ginv(w`*w); trouble2=w`*y;
				resids1to2[, j]=y-w*trouble1*trouble2; /* Use X2 but B2-tilde (obtained by estimating B2 using X2 while assuming B1-hat structure) to find residuals */
				err21[, j]=y-w*bigBnow1[whichcol, j];  /* Use X2 data but B1-hat (obtainedestimating B1 using X1) to find residuals */

				yw=allx1[whichrow, j]||allx1[whichrow, whichcol];
				yw=yw-yw[:,]; 
				y=yw[, 1]; w=yw[, 2:(d+1)];
/*				trouble1=ginv(w`*w); trouble2=w`*y;*/
/*				resids1to1[, j]=y-w*trouble1*trouble2;*/
				resids1to1[, j]=y-w*bigBnow1[whichcol, j];
			end;
			else do; 
				y=allx2[whichrow, j];
				resids1to2[, j]=y-y[:];
				err21[, j]=y-y[:]; 

				y=allx1[whichrow, j];
				resids1to1[, j]=y-y[:];
			end;
		
			/* second assume G2 structure */
			whichcol=loc(bigBnow2[, j]^=0); 

			d=ncol(whichcol);
			if (d>0) then do;
				/*~~~~~~~~~~~~ center but not standardize w and y ~~~~~~~~~~~~*/
				yw=allx1[whichrow, j]||allx1[whichrow, whichcol];
				yw=yw-yw[:,]; 
				y=yw[, 1]; w=yw[, 2:(d+1)];
				trouble1=ginv(w`*w); trouble2=w`*y;
				resids2to1[, j]=y-w*trouble1*trouble2;
				err12[, j]=y-w*bigBnow2[whichcol, j]; 

				yw=allx2[whichrow, j]||allx2[whichrow, whichcol];
				yw=yw-yw[:,]; 
				y=yw[, 1]; w=yw[, 2:(d+1)];
/*				trouble1=ginv(w`*w); trouble2=w`*y;*/
/*				resids2to2[, j]=y-w*trouble1*trouble2;*/
				resids2to2[, j]=y-w*bigBnow2[whichcol, j]; 
			end;
			else do; 
				y=allx1[whichrow, j];
				resids2to1[, j]=y-y[:];
				err12[, j]=y-y[:]; 

				y=allx2[whichrow, j];
				resids2to2[, j]=y-y[:];
			end;
		end;

		firsthalf=floor(dall/2); 
		samsize1=nj*firsthalf; samsize2=NoObs-samsize1;
		pvs1to2=j(NoNode, 3, 1); pvs2to1=j(NoNode, 3, 1); 
		pvs1to1=j(NoNode, 3, 1); pvs2to2=j(NoNode, 3, 1); 
		disc=j(NoNode, 1, 1);
		PIscore=0; tophowmany=0; 
		do sp=1 to split;
			mark=j(NoObs, 1, 2);
			call randseed(randomseed+sp);
			twoset=ranperm(dall);
			do expcond=1 to firsthalf;
				thiscond=twoset[expcond];
				mark[((thiscond-1)*nj+1):(thiscond*nj)]=1; 
			end;
			resids1=resids1to2[loc(mark=1), ];
			resids2=resids1to2[loc(mark=2), ];
			residsrv1=resids2to1[loc(mark=1), ];
			residsrv2=resids2to1[loc(mark=2), ];

			g1resids1=resids1to1[loc(mark=1), ]; 
			g1resids2=resids1to1[loc(mark=2), ]; 
			g2resids1=resids2to2[loc(mark=1), ]; 
			g2resids2=resids2to2[loc(mark=2), ]; 

			do j=1 to NoNode;
				/* Going from G1 to G2 */
				sample1=resids1[, j];
				sample2=resids2[, j];
				howcomp=meanvartest(sample1, sample2, samsize1, samsize2);
				pvs1to2[j, ]=howcomp;

				/* Going from G2 to G1 */
				sample1=residsrv1[, j];
				sample2=residsrv2[, j];
				howcomp=meanvartest(sample1, sample2, samsize1, samsize2);  
				pvs2to1[j, ]=howcomp;

				/* Going from G1 to G1 */
				sample1=g1resids1[, j];
				sample2=g1resids2[, j];
				howcomp=meanvartest(sample1, sample2, samsize1, samsize2);
				pvs1to1[j, ]=howcomp;

				/* Going from G2 to G2 */
				sample1=g2resids1[, j];
				sample2=g2resids2[, j];
				howcomp=meanvartest(sample1, sample2, samsize1, samsize2);
				pvs2to2[j, ]=howcomp;


				/* Computer DISCERN score */
				if sp=1 then do;
					commondeno=resids1to1[,j]`*resids1to1[,j]+resids2to2[,j]`*resids2to2[,j];
					disc[j]=(resids2to1[,j]`*resids2to1[,j]+resids1to2[,j]`*resids1to2[,j])/commondeno;
				end; 
			end;

			PIscore=PIscore+(exp(-pvs1to2[, 3])<>exp(-pvs2to1[, 3]))/(exp(-pvs1to1[, 3])<>exp(-pvs2to2[, 3]));
			tophowmany=tophowmany+(sum(pvs1to2[, 3]<pvcutoff))<>(sum(pvs2to1[, 3]<pvcutoff)); /* This gives better results than below */
/*			tophowmany=tophowmany+sum((pvs1to2[, 3]><pvs2to1[, 3])<pvcutoff);*/
		end;
		PIscore=PIscore/split;
		tophowmany=ceil(tophowmany/split);
finish; 


/* Benjamini-Hochberg procedure: input original pvalues, output the node indices with significant pvalues after BH adjustment */
start BHpv(pvs) global(fdr); 
	m=max(nrow(pvs),ncol(pvs)); 
	ids=do(1, m, 1)`;
	tag=ids||pvs;
	call sort(tag, 2);
	crt=fdr*(ids/m); /* this is the BH critical value for each sorted pvalue */
	lower=loc(tag[, 2]<crt);
	if ncol(lower)>0 then do;
		cutoff=lower[ncol(lower)];
		claim=tag[, 1][1:cutoff];
		return(claim); 
	end;
	else do;
		return(-1); /* this means no significant pvalues */
	end;
finish; 

/* Find driver nodes after identifying differential node */
start FindDriver(DiffNode, bigBnow1, bigBnow2) global(drivernodes, nodriv, notdriv);
	tophowmany=max(nrow(DiffNode), ncol(DiffNode)); 
	NoNode=ncol(bigBnow1); 
	driver=j(tophowmany, NoNode, 0);
	do k=1 to tophowmany;
		child=DiffNode[k]; 
		pa1=loc(bigBnow1[, child]^=0);
		pa2=loc(bigBnow2[, child]^=0);
		nopa1=ncol(pa1); nopa2=ncol(pa2); 

		if ((nopa1><nopa2)>0) then do;
			discrep=setdif(union(pa1, pa2), xsect(pa1, pa2)); /* the symmetric difference between two parent sets */
			if ncol(discrep)>0 then driver[k, discrep]=1;
		end;
		else do; 
			if (nopa1=0 & nopa2>0) then driver[k, pa2]=1; 
			if (nopa2=0 & nopa1>0) then driver[k, pa1]=1; 
		end;
	end;
	touched=loc(driver[+, ]>0);
	claimdriv=ncol(touched);
/*	if claimdriv>0 then do;*/
		pickdriv=sum(element(drivernodes, touched)); 
		tpr4driv=pickdriv/nodriv; 
		spe4driv=ncol(xsect(setdif(1:NoNode, touched), notdriv))/(NoNode-nodriv);  
		fdr4driv=ncol(xsect(touched, notdriv))/claimdriv; 
/*	end;*/
*	else do;  /* I do not claim any driver node */
/*		tpr4driv=0; spe4driv=1; fdr4driv=0; touched=j(1, NoNode, -88); */
/*	end;*/
	return(tpr4driv||spe4driv||fdr4driv||touched);
finish;

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
pall=30; rho1=0.3; rho2=0.3; bb=pall*rho1; cc=pall*rho2; aa=pall-bb-cc; 
nodiff=bb; wanttotal=4*pall;
ids=do(1, pall, 1)`; 
dall=pall-1; diagIdx = do(1, pall*pall, pall+1);	
sige=sqrt(0.25); sigx=1; 
optn=(pall-1)||0;
allpossible=pall*(pall-1)/2; 
allneg=allpossible-wanttotal;

pi=constant("pi"); twopi=2*pi; 
/*sndelta=sqrt(pi*(1-sige#sige)/2); */
sndelta=0.9;
mnshift=sqrt(2/pi)*sndelta;
NoNode=pall; 
dall=NoNode-1;

maxit=50; small=1e-3; small2=1e-4; small3=1e-6; threshold=1e-4; eps=1e-20; 
scada=3.7; 
varpvct=0.05; fdr=0.05;
pvcutoff=0.025;

lambmin=0.01; lambmax=0.5;
pvout=j(1, pall, .); 

rt=lambmin/lambmax;
m=20;
q=rt**(1/(m-1));

permno=300; 
/*njs={5, 10, 20, 30, 40, 50, 60}; */
njs={5, 10}; 
nonj=nrow(njs);
topds=do(1, dall, 1); 

testout=j(1, 3+2*2*2, 0);

nographs=1;
MCrep=100; jumpstart=4;
split=10; case=1;
if case=1 then do; 
	maxparent=6; wanttotal=4#pall; 
end;
 
if case=2 then do;
	maxparent=20;
end;
takethisseed=0; 
do nog=1 to nographs;
	/* create two regression coefficients matrices */
	noedges=0; 
	edseed=takethisseed+1;
	blank=1; 

	/*~~~~~~~~~~~~~ Case 1: aa non-diff/non-driver nodes, bb diff nodes, and at most cc driver nodes ~~~~~~~~~~~~~*/
	if case=1 then do;
		do while (noedges^=wanttotal | blank=1);
			edseed=edseed+1; 
			call randseed(edseed);
			whatorder=ranperm(pall);
			bigBx=j(pall, pall, 0); 
			bigBx[, whatorder[1]]=0;
			nopar=j(pall, 1, 0);
			do k=2 to pall;
				thisnode=whatorder[k]; 
				nopar[k]=sample(0:min((k-1), maxparent), 1);
				if nopar[k]>0 then do;
					whichpar=sample(whatorder[1:k-1], nopar[k], "NoReplace"); 
					bigBx[whichpar, thisnode]=1;
				end; 
			end;
			noedges=ncol(loc(bigBx^=0)); 
			if noedges=wanttotal then do; 
				block3=bigBx[(aa+bb+1):pall, (aa+1):(aa+bb)];
				if (min(block3[+, ])>0) then blank=0; 
			end;
		end;
		takethisseed=edseed;
		
		/* Create another graph by deleting nonzero entries in certain blocks of bigBx  */
		bigBy=bigBx; 
		block3=bigBy[(aa+bb+1):pall, (aa+1):(aa+bb)];
		drivernodes=0;
		do k=1 to bb;
			candidates=loc(block3[, k]^=0); 
			call randseed(k+edseed);
			driver=sample(candidates, 1);
			block3[driver, k]=0;
			drivernodes=drivernodes||(driver+aa+bb);
		end;
		bigBy[(aa+bb+1):pall, (aa+1):(aa+bb)]=block3;
		noedgesy=sum(bigBy^=0);
		print noedges noedgesy;
		
		diffnodes=do(aa+1, aa+bb, 1); 
		drivernodes=unique(drivernodes[2:ncol(drivernodes)]);

		diffonly=diffnodes; 
		driveronly=drivernodes; 
		
		notdiff=setdif(ids, diffnodes); 
		notdriv=setdif(ids, drivernodes);
 
		nodriv=ncol(drivernodes);
		nodiff=ncol(diffnodes);
		print diffnodes drivernodes;	
	end;

	
	/*~~~~~~~~~~~~~~~ Case 2: at most bb diff nodes, at most cc driver nodes, allow overlapping ~~~~~~~~~~~~~~~~*/
	if case=2 then do;
		blank13=1;
		do while (blank13=1);
			edseed=edseed+1; 
			call randseed(edseed);
			whatorder=ranperm(pall);
			bigBx=j(pall, pall, 0); 
			bigBx[, whatorder[1]]=0;
			nopar=j(pall, 1, 0);
			do k=2 to pall;
				thisnode=whatorder[k]; 
				nopar[k]=sample(0:min((k-1), maxparent), 1);
				if nopar[k]>0 then do;
					whichpar=sample(whatorder[1:k-1], nopar[k], "NoReplace"); 
					bigBx[whichpar, thisnode]=1;
				end; 
			end;
			noedges=ncol(loc(bigBx^=0));

			block13=bigBx[(aa+1):pall, (aa+1):(aa+bb)];
			if (min(block13[+, ])>0) then blank13=0;
		end;
		takethisseed=edseed;
		wanttotal1=noedges; 

		/* Create another graph by deleting nonzero entries in certain blocks of bigBx  */
		bigBy=bigBx; 
		block3=bigBy[(aa+bb+1):pall, (aa+1):(aa+bb)];
		block1=bigBy[(aa+1):(aa+bb), (aa+1):(aa+bb)];
		drivernodes=0; diffnodes=0;
		do k=1 to bb;
			candidates=loc(block3[, k]^=0);
			if ncol(candidates)>0 then do;
				call randseed(k);
				driver=sample(candidates, 1);
				block3[driver, k]=0;
				diffnodes=diffnodes||(k+aa);
				drivernodes=drivernodes||(driver+aa+bb);
			end; 

			candidates=loc(block1[, k]^=0); 
			if ncol(candidates)>0 then do;
				call randseed(k+9000);
				driver=sample(candidates, 1);
				block1[driver, k]=0;
				diffnodes=diffnodes||(k+aa);
				drivernodes=drivernodes||(driver+aa);
			end;
		end;
		bigBy[(aa+bb+1):pall, (aa+1):(aa+bb)]=block3;
		bigBy[(aa+1):(aa+bb), (aa+1):(aa+bb)]=block1;
		noedges2=ncol(loc(bigBy^=0));
		wanttotal2=noedges2;

		print noedges noedges2;
		diffnodes=unique(diffnodes[2:ncol(diffnodes)]);
		drivernodes=unique(drivernodes[2:ncol(drivernodes)]);
		bothnodes=xsect(diffnodes, drivernodes);
		driveronly=setdif(drivernodes, bothnodes);

		diffonly=setdif(diffnodes, bothnodes); 
		notdiff=setdif(ids, diffnodes); 
		notdriv=setdif(ids, drivernodes);
		
		nodriv=ncol(drivernodes);
		nodiff=ncol(diffnodes); 
		noboth=ncol(bothnodes);
		print diffnodes drivernodes bothnodes;
	end;
	neitherrole=setdif(setdif(1:pall, diffnodes), drivernodes); 

	u=j(pall*pall, 1, 1);
	call randseed(nog+77777);
	call randgen(u, "Uniform"); 
	
	fillin=loc(bigBx^=0);
	call randseed(nog+333321); 
	bf = rand("Bernoulli", repeat(0.5, ncol(fillin)));
	bigBx[fillin]=u[1:ncol(fillin)]+bf-(1-bf)*2;

	call randseed(nog+765);
	call randgen(u, "Uniform"); 
	
	fillin=loc(bigBy^=0);
	call randseed(nog+321); 
	bf = rand("Bernoulli", repeat(0.5, ncol(fillin)));
	bigBy[fillin]=u[1:ncol(fillin)]+bf-(1-bf)*2;

	revisey=Khan(bigBy);
	bigBy=revisey[, 1:pall]; 
	whatordery=revisey[, pall+1]`;
	nopary=j(pall, 1, 0);
	do k=2 to pall; 
		thischild=whatordery[k];
		itsparent=loc(bigBy[, thischild]^=0);
		if ncol(itsparent)>0 then nopary[k]=ncol(itsparent);
	end; 

	
	do whichnj=1 to nonj; 
		nj=njs[whichnj]; totn=pall*nj;  
		NoObsInt=totn; NoObs=totn-nj;
		bicpen1=log(NoObs)/NoObs; 
		NoObsInt2=NoObsInt; NoObs2=totn-nj;
		bicpen2=log(NoObs2)/NoObs2;
		
		do rep=1 to MCrep;  
			duetoemp:
		   	seedx=j(NoObs, 1, 341+nog+rep+whichnj+jumpstart); 
			seedxj=j(nj, 1, 451+nog+rep+whichnj+jumpstart);
			seede=j(NoObs, 1, 4351+nog+rep+whichnj+jumpstart); 
			seedu=j(NoObsInt, pall, 321+nog+rep+whichnj+jumpstart);	
			
			allx1=j(totn, pall, 0); allx2=allx1; 
			call randseed(399999+nog+rep+whichnj+jumpstart);
			u = j(2*nj, 1);                
			call randgen(u, "Uniform"); 
			
			root=whatorder[1]; 
			if (root>1 & root<pall) then whichrow=(1:((root-1)#nj))||((root#nj+1):totn);
			if root=1 then whichrow=(nj+1):totn;
			if root=pall then whichrow=1:((pall-1)#nj);
			allx1[((root-1)*nj+1):(root*nj), whatorder[1]]=u[1:nj]; /* interventional data for the root node */
			allx1[whichrow, whatorder[1]]=sige#normal(seedx);    /* not interventional data but data according to the regression model for a root node */

			root=whatordery[1]; 
			if (root>1 & root<pall) then whichrow=(1:((root-1)#nj))||((root#nj+1):totn);
			if root=1 then whichrow=(nj+1):totn;
			if root=pall then whichrow=1:((pall-1)#nj);
			allx2[((root-1)*nj+1):(root*nj), whatordery[1]]=u[(nj+1):(2*nj)]; /* interventional data for the root node */
			allx2[whichrow, whatordery[1]]=sige#normal(seedx+9754); /* not interventional data but data according to the regression model for a root node */
			do k=2 to pall;
				thisnode=whatorder[k]; 
				call randseed(123+nog+rep+whichnj+k);
				u = j(2*nj, 1);
				call randgen(u, "Uniform"); 
				if (thisnode>1 & thisnode<pall) then whichrow=(1:((thisnode-1)#nj))||((thisnode#nj+1):totn);
				if thisnode=1 then whichrow=(nj+1):totn;
				if thisnode=pall then whichrow=1:((pall-1)#nj);
					
				allx1[((thisnode-1)#nj+1):(thisnode#nj), thisnode]=u[1:nj];
				if nopar[k]>0 then do;
					whichpar=loc(bigBx[, thisnode]^=0);
					/* the following is to generate Gaussian model error */
					allx1[whichrow, thisnode]=allx1[whichrow, whichpar]*bigBx[whichpar, thisnode]+sige#normal(seede+k);

					/* the following is to generate skewed normal model error */
/*					modelerr1=normal(seede+k); modelerr2=normal(seede+k+99999);*/
/*					modelerr=sndelta#modelerr1+sqrt(1-sndelta#sndelta)#modelerr2;*/
/*					modelerr=modelerr#(modelerr1>=0)-modelerr#(modelerr1<0);*/
/*					modelerr=modelerr-mnshift;*/
/*					allx1[whichrow, thisnode]=sqrt(abs(allx1[whichrow, whichpar]))*bigBx[whichpar, thisnode]+modelerr;*/
				end; 
				else allx1[whichrow, thisnode]=sige#normal(seedx+k);
				
				thisnode=whatordery[k]; 
				if (thisnode>1 & thisnode<pall) then whichrow=(1:((thisnode-1)#nj))||((thisnode#nj+1):totn);
				if thisnode=1 then whichrow=(nj+1):totn;
				if thisnode=pall then whichrow=1:((pall-1)#nj);
					
				allx2[((thisnode-1)#nj+1):(thisnode#nj), thisnode]=u[(nj+1):(2*nj)];
				if nopary[k]>0 then do;
					whichpar=loc(bigBy[, thisnode]^=0);
					/* generate Gaussian DAG */
/*					if (thisnode=15) then */
/*					allx2[whichrow, thisnode]=sqrt(abs(allx2[whichrow, whichpar]))*bigBy[whichpar, thisnode]+sige#normal(seede+k+8766655);*/
/*					else */
					allx2[whichrow, thisnode]=allx2[whichrow, whichpar]*bigBy[whichpar, thisnode]+sige#normal(seede+k+8766655);
					/* the following is to generate skewed normal model error */
/*					modelerr1=normal(seede+k+6666666); modelerr2=normal(seede+k+33333);*/
/*					modelerr=sndelta#modelerr1+sqrt(1-sndelta#sndelta)#modelerr2;*/
/*					modelerr=modelerr#(modelerr1>=0)-modelerr#(modelerr1<0);*/
/*					modelerr=modelerr-mnshift;*/
/*					allx2[whichrow, thisnode]=allx2[whichrow, whichpar]*bigBy[whichpar, thisnode]+modelerr;*/
				end; 
				else allx2[whichrow, thisnode]=sige#normal(seedx+k+9898765);
			end;

			NoObs=totn-nj;
	
			/*---------------- Done with graph and data generation. ----------------*/
			bigBnow1=j(pall, pall, 0);
			/**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Estimate B1 & B2 ~~~~~~~~~~~~~~~~~~~~~~~~~~**/
			/* Get the non-penalized ls score estimates as starting values of bigB. */
			bigBini1=j(NoNode, NoNode, 0); bigBini2=j(NoNode, NoNode, 0); 
			do j=1 to NoNode; 
				ignorej=loc((1:NoNode)^=j);
				if (j>1 & j<NoNode) then whichrow=(1:((j-1)#nj))||((j#nj+1):NoObsInt);
				if j=1 then whichrow=(nj+1):NoObsInt;
				if j=NoNode then whichrow=1:((NoNode-1)#nj);

				w=allx1[whichrow, ignorej];
				y=allx1[whichrow, j];
				y=y-y[:]; w=w-w[:, ]; 
				bigBini1[ignorej, j]=inv(w`*w)*w`*y;

				w=allx2[whichrow, ignorej];
				y=allx2[whichrow, j];
				y=y-y[:]; w=w-w[:, ]; 
				bigBini2[ignorej, j]=inv(w`*w)*w`*y;
			end;
			
			/* make the initial B not contain [i,j] and [j,i] entries both nonzero */
			do i=1 to (NoNode-1); 
				do j=i+1 to NoNode;
					if abs(bigBini1[i,j])<=abs(bigBini1[j,i]) then bigBini1[i,j]=0;
					else bigBini1[j,i]=0;

					if abs(bigBini2[i,j])<=abs(bigBini2[j,i]) then bigBini2[i,j]=0;
					else bigBini2[j,i]=0;
				end;
			end;

			run Tuning(lambda1, lambda2, bigBnow1, bigBnow2, lambmax, m, q, allx1, allx2, bigBini1, bigBini2);

			/******************************* Use tuning parameter selected above to get graphs *******************************/
			randomseed=34646+nog+rep+whichnj;
			run Find2score(PIscore, disc, tophowmany, allx1, allx2, bigBnow1, bigBnow2, randomseed); 
	
			if (norm(bigBnow1)=0 | norm(bigBnow2)=0) then do;
				jumpstart=jumpstart+77777;
				goto duetoemp;
			end;

			sortedPI=ids||PIscore;
			call sort(sortedPI, 2, 2); 
			sortedDISC=ids||disc;
			call sort(sortedDISC, 2, 2); 
					
			do whichtopd=1 to dall;
				dop=topds[whichtopd];
				/* Claim the top "dop" nodes as differential nodes based on the score ranking */	
				diffnd=sortedPI[1:dop, 1];
				pickdiff=sum(element(diffnodes, diffnd));  /* the number of differential nodes that are also the top tophowmany nodes according to PI score */
				tpr4diff=pickdiff/nodiff; 
				fpr4diff=ncol(xsect(diffnd, notdiff))/(NoNode-nodiff); 
				PIdiff_rnk=tpr4diff||fpr4diff; 
				allout=FindDriver(diffnd, bigBnow1, bigBnow2);
				PIdriv_rnk=allout[1]||(1-allout[2]); /* record tpr and fpr for driver nodes identification */
				
				diffnd=sortedDISC[1:dop, 1];
				pickdiff=sum(element(diffnodes, diffnd));  /* the number of differential nodes that are also the top tophowmany nodes according to DISCERN2 score */
				tpr4diff=pickdiff/nodiff; 							
				fpr4diff=ncol(xsect(diffnd, notdiff))/(NoNode-nodiff); 
				discdiff_rnk=tpr4diff||fpr4diff;
				allout=FindDriver(diffnd, bigBnow1, bigBnow2);
				discdriv_rnk=allout[1]||(1-allout[2]);
				testout=testout//(nj||tophowmany||dop||PIdiff_rnk||discdiff_rnk||PIdriv_rnk||discdriv_rnk);			
			end;
		end;
	end;
end;
testout=testout[2:nrow(testout), ];

create out.roc08082023c1 from testout[colname=('nj'||'tophowmany'||'topd'
||'PItpr4diff_rnk'||'PIfpr4diff_rnk'||'DISCtpr4diff_rnk'||'DISCfpr4diff_rnk'
||'PItpr4driv_rnk'||'PIfpr4driv_rnk'||'DISCtpr4driv_rnk'||'DISCfpr4driv_rnk')];
append from testout; 

quit; 
proc sort data=out.roc08082023c1; by nj topd;
run; 

proc means data=out.roc08082023c1 mean stderr median n;
	by nj topd;
run;


data track;
    t=time(); d=today();
    put 'Time: ' t time9.;
    put 'Date: ' d date9.;
run;

quit;

