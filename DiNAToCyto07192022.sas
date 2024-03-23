dm "log; clear;";
options ps=80 ls=90 pageno=1 nodate formdlim='-' notes threads;
title;
ODS HTML CLOSE;
ODS HTML;
/*************************************************************************************/
/* Date: 07-19-2022                                                                  */
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
/*       (107) On 05-25-2021 revised DiNA05032021.sas to apply to flow cytometry data*/
/*       (108) On 01-26-2022 revised DiNAToCyto05252021.sas and DiNA10112021.sas to  */
/*             apply to flow cytometry data.                                         */
/*       (109) On 05-05-2022 revised DiNAToCyto01262022.sas to add competing methods */
/*             based on DISCERN idea.                                                */
/*       (110) On 07-19-2022 revised DiNAToCyto05052022.sas to add estimating B based*/
/*             on permuted data.                                                     */
/*************************************************************************************/
/*libname out "Z:\Research_Networks\Output";*/
/*libname out "C:\Research_Networks\Output";*/
/*libname out "E:\Research_Networks";*/
/*libname out "C:\Documents and Settings\huang\Research_Networks";*/
libname rawdata "C:\Research_Networks\Data";

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


/* Use my earlier variable selection code to select parents for each node based on corrected score method with LQA algorithm*/
start findgraph2(lambda, allw, bigBini0, mark) global(dall, pall, threshold, maxit, small, small2, scada);
	bigBnow=bigBini0; itno=0; absdiff=10; 
	n=nrow(allw);
	do while (absdiff>threshold & itno<=maxit);
		bigBlast=bigBnow; itno=itno+1;
 
		do j=1 to pall;
			ignorej=loc((1:pall)^=j);
			betanow=bigBnow[, j];
			
			if max(abs(betanow))=0 then do;
				rc=-1; goto doneold; 
			end; 

			whichrow=loc(mark[, j]=1);
			n=sum(mark[,j]);  
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
	return(bigBnow); 
finish;


/* obtain a p-value matrix for the unpenalized estimated regression coefficients obtained from corrected score method */
start getpv2(graphest, allw, mark) global(n, y, w, d);
	pall=ncol(allw); totn=nrow(allw);  
	pvmat=j(pall, pall, -8); 
	do j=1 to pall;
		whichcol=loc(graphest[, j]^=0); d=ncol(whichcol); 
		if (d>0) then do;
			whichrow=loc(mark[, j]=1);
			nall=sum(mark[,j]);
		
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
	return(graphest);
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

/* Find driver nodes after identifying node */
start FindDriver(DiffNode) global(pall, bigBnow1, bigBnow2);
	tophowmany=max(nrow(DiffNode), ncol(DiffNode)); 
	driver=j(tophowmany, pall, 0);
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
	if claimdriv>0 then return(touched);
	else return(-1); 
finish;


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
use rawdata.all9; read all into alldata; close rawdata.all9;
pall=11;

/*~~~~~~~~~~~~~~~~~~ if I want to pick out only nodes that have both observational data and interventional data, do the following ~~~~~~~~~~~~~~~~~~*/
/*goodnodes=loc((alldata[, 1:11])[+, ]<nrow(alldata));*/
goodnodes=1:pall;
print goodnodes;

alldata=alldata[, goodnodes||(goodnodes+11)];
pall=ncol(goodnodes); 
/*print "No. of nodes with both observational data and interventional data:" pall; */
/*~~~~~~~~~~~~~~~~~~ end of picking out nodes ~~~~~~~~~~~~~~~~~~*/

dall=pall-1; 
optn=dall||0;

/* There are 9 experimental conditions in this data. */
intervene_size={853, 902, 911, 723, 810, 799, 848, 913, 707};
exp_cond=nrow(intervene_size); 
stacknj=cusum(intervene_size);

/*----------------------------------------------------------------------------------------------------*/
/*~~~~~~~~~~~~~~~~~~ if I only use part of the full data, do the following ~~~~~~~~~~~~~~~~~~*/
reduct=1; /*set reduct bigger than 1 if I want to keep halving the full data */
do re=1 to reduct; 
	newintervene_size=round(intervene_size/2);
	newstacknj=cusum(newintervene_size);
	newalldata=alldata[1:newintervene_size[+], ];
	do itv=1 to exp_cond; 
		if itv=1 then do;
			call randseed(itv+re+987876);
			permindex=ranperm(intervene_size[itv]);

			messed=alldata[permindex, ];
			newalldata[1:newintervene_size[1], ]=messed[1:newintervene_size[1], ];
		end;	
		else do;
			call randseed(itv+re+677);
			permindex=ranperm(intervene_size[itv])+stacknj[itv-1];

			messed=alldata[permindex, ];
			newalldata[(newstacknj[itv-1]+1):newstacknj[itv], ]=messed[1:newintervene_size[itv], ];
		end;
	end;
	alldata=newalldata; 
	intervene_size=newintervene_size;
	exp_cond=nrow(intervene_size); 
	stacknj=cusum(intervene_size);
end;
/*~~~~~~~~~~~~~~~~~~ end of picking out a subset of the full data ~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*----------------------------------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*~~~~~~~~~~~~~~~~~~ if I split the full data into two halves, do the following ~~~~~~~~~~~~~~~~~~*/
/* Later split the whole data into two data sets. The following record the number of data in each of the nine experimental conditions for each data set*/
/*intervene_cond1=round(intervene_size/2);*/
/*intervene_cond2=intervene_size-intervene_cond1;*/
/*stacknj1=cusum(intervene_cond1);*/
/*stacknj2=cusum(intervene_cond2);*/
/*alldata1=j(stacknj1[exp_cond], 2*pall, 0); */
/*alldata2=j(stacknj2[exp_cond], 2*pall, 0); */
/**/
/*do itv=1 to exp_cond; */
/*	if itv=1 then do;*/
/*		cond1_idx=j(intervene_cond1[1], 1, 1);*/
/*		cond2_idx=j(intervene_cond2[1], 1, 1);*/
/*		*/
/*		call randseed(itv);*/
/*		permindex=ranperm(intervene_size[itv]);*/
/**/
/*		messed=alldata[permindex, ];*/
/*		alldata1[1:intervene_cond1[1], ]=messed[1:stacknj1[1], ];*/
/*		alldata2[1:intervene_cond2[1], ]=messed[(stacknj1[1]+1):stacknj[1], ];	*/
/*	end;	*/
/*	else do;*/
/*		cond1_idx=cond1_idx//j(intervene_cond1[itv], 1, itv);*/
/*		cond2_idx=cond2_idx//j(intervene_cond2[itv], 1, itv);*/
/*		 */
/*		call randseed(itv);*/
/*		permindex=ranperm(intervene_size[itv])+stacknj[itv-1];*/
/**/
/*		messed=alldata[permindex, ];*/
/*		alldata1[(stacknj1[itv-1]+1):stacknj1[itv], ]=messed[1:intervene_cond1[itv], ];*/
/*		alldata2[(stacknj2[itv-1]+1):stacknj2[itv], ]=messed[(intervene_cond1[itv]+1):intervene_size[itv], ];*/
/*	end;*/
/*end;*/

/* The following alldata results from permuting each of the nine excel-file data set (under nine experimental conditions), then combine the nine permuted data */
/*mark1=alldata1[, 1:pall];*/
/*allx1=alldata1[, (pall+1):2*pall];*/
/*mark2=alldata2[, 1:pall];*/
/*allx2=alldata2[, (pall+1):2*pall];*/
/*~~~~~~~~~~~~~~~~~~ end of creating two data sets by splitting the full data ~~~~~~~~~~~~~~~~~~*/
/*---------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*~~~~~~~~~~~~~~~~~~ if I want to keep the original (part of) full data and then create another data of similar structure by manipulating the full data, do the following ~~~~~~~~~~~~~~~~~~*/
mark1=alldata[, 1:pall];
allx1=alldata[, (pall+1):2*pall]; 
mark2=mark1;
noise=2*normal(j(nrow(allx1), ncol(allx1), 4656));
allx2=allx1+noise; 


/* standardize each column */
eachmean=mean(allx1); eachstd=std(allx1);
allx1=(allx1-eachmean@j(nrow(allx1), 1, 1))/(eachstd@j(nrow(allx1), 1, 1));
eachmean=mean(allx2); eachstd=std(allx2);
allx2=(allx2-eachmean@j(nrow(allx2), 1, 1))/(eachstd@j(nrow(allx2), 1, 1));


messobs1=7; /* the node whose observational data I will mess up in the second data set (in order to make up a differential node) */
*messobs2=7; /* the node whose observational data I will mess up in the second data set */
messint=2; /* the node whose interventional data I will mess up in the second data set (maybe to create a driver node) */
allx2[loc(mark2[, messobs1]=1), messobs1]=median(allx2[loc(mark2[, messobs1]=1), messobs1])+normal(j(ncol(loc(mark2[, messobs1]=1)), 1, 5656));
/*allx2[loc(mark2[, messobs2]=1), messobs2]=median(allx2[loc(mark2[, messobs2]=1), messobs2])+normal(j(ncol(loc(mark2[, messobs2]=1)), 1, 447));*/
/*allx2[loc(mark2[, messint]=0), messint]=median(allx2[loc(mark2[, messint]=0), messint])+normal(j(ncol(loc(mark2[, messint]=0)), 1, 656));*/
messintsch=8;
allx2[loc(mark2[, messint]=0), messint]=median(allx2[loc(mark2[, messintsch]=1), messintsch])+normal(j(ncol(loc(mark2[, messint]=0)), 1, 656));
/*allx2[loc(mark2[, messint]=0), messint]=median(allx2[loc(mark2[, messint]=0), messint]);*/

intervene_cond1=intervene_size;
intervene_cond2=intervene_cond1;
stacknj1=cusum(intervene_cond1);
stacknj2=cusum(intervene_cond2);

do itv=1 to exp_cond; 
	if itv=1 then do;
		cond1_idx=j(intervene_cond1[1], 1, 1);
		cond2_idx=j(intervene_cond2[1], 1, 1);
	end;	
	else do;
		cond1_idx=cond1_idx//j(intervene_cond1[itv], 1, itv);
		cond2_idx=cond2_idx//j(intervene_cond2[itv], 1, itv);
	end;
end;
/*~~~~~~~~~~~~~~~~~~ end of creating another data set by manipulating the original full data ~~~~~~~~~~~~~~~~~~*/
/*----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/



sampsize1=nrow(allx1); sampsize2=nrow(allx2);
print sampsize1 sampsize2;

dall=pall-1; diagIdx = do(1, pall*pall, pall+1);	
optn=(pall-1)||0;
pi=constant("pi"); twopi=2*pi;

maxit=100; small=1e-3; small2=1e-4; small3=1e-6; threshold=1e-4; eps=1e-20; 
scada=3.7; 
varpvct=0.05; fdr=0.05;

pvcutoff=0.025;

rt=0.2;
lambmax=1;
lambmin=lambmax*rt;
m=100; 
aa=0.1; scada=3.7; scada1=scada-1;
q=rt**(1/(m-1));
split=2; 
permno=300; 

/**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Estimate B1 & B2 ~~~~~~~~~~~~~~~~~~~~~~~~~~**/
/* Get the non-penalized ls score estimates as starting values of bigB. */
NoNode=ncol(allx1); NoObsInt=nrow(allx1); 
NoObs1=j(NoNode, 1, 0); NoObs2=j(NoNode, 1, 0); 
obssize1=j(NoNode, 1, 0); obssize2=j(NoNode, 1, 0);
dall=NoNode-1; 
bigBini1=j(NoNode, NoNode, 0); bigBini2=j(NoNode, NoNode, 0); 
do j=1 to NoNode; 
	ignorej=loc((1:NoNode)^=j);
	whichrow=loc(mark1[, j]=1);
	NoObs1[j]=sum(mark1[,j]);  
	obssize1[j]=sum(mark1[, j]); obssize2[j]=sum(mark2[, j]);
			
	w=allx1[whichrow, ignorej];
	y=allx1[whichrow, j];
	y=y-y[:]; w=w-w[:, ]; 
	bigBini1[ignorej, j]=inv(w`*w)*w`*y;

	whichrow=loc(mark2[, j]=1);
	NoObs2[j]=sum(mark2[,j]);
	w=allx2[whichrow, ignorej];
	y=allx2[whichrow, j];
	y=y-y[:]; w=w-w[:, ]; 
	bigBini2[ignorej, j]=inv(w`*w)*w`*y;
end;

bicpen=(log(NoObs1)/NoObs1)[:];
bicpen2=(log(NoObs2)/NoObs2)[:];

/* make the initial B not contain [i,j] and [j,i] entries both nonzero */
do i=1 to (NoNode-1); 
	do j=i+1 to NoNode;
		if abs(bigBini1[i,j])<=abs(bigBini1[j,i]) then bigBini1[i,j]=0;
		else bigBini1[j,i]=0;

		if abs(bigBini2[i,j])<=abs(bigBini2[j,i]) then bigBini2[i,j]=0;
		else bigBini2[j,i]=0;
	end;
end;

/* Estimate bigB for a sequence of candidate lambda's */  
lambda1=lambmax; lambda2=lambmax;
do mm=1 to m;
	lambda=lambmax*q**(mm-1); 
	
	bigBnow1=findgraph2(lambda, allx1, bigBini1, mark1); 
	pvmat=getpv2(bigBnow1, allx1, mark1);
	bigBnow1=topsort(bigBnow1, pvmat);
	
	bigBnow2=findgraph2(lambda, allx2, bigBini2, mark2); 
	pvmat=getpv2(bigBnow2, allx2, mark2);
	bigBnow2=topsort(bigBnow2, pvmat);
	
	/* using the estimated graph, compute SIC tuning parameter selector (eveluated at the penalized Bhat). */
	sic1=0; sic2=0;
	do j=1 to NoNode; 
		ignorej=loc((1:NoNode)^=j); 
		whichrow=loc(mark1[, j]=1);
		
		yw=allx1[whichrow, ]; 
		yw=yw-yw[:,];
		y=yw[, j]; w=yw[, ignorej];
		sic1=sic1+lsobj(bigBnow1[ignorej, j], y, w);

		whichrow=loc(mark2[, j]=1);
		
		yw=allx2[whichrow, ]; 
		yw=yw-yw[:,];
		y=yw[, j]; w=yw[, ignorej];
		sic2=sic2+lsobj(bigBnow2[ignorej, j], y, w);
	end;

	if mm=1 then do; 
		smallestsic1=sic1+sum(bigBnow1^=0)#bicpen; 
		sic1=smallestsic1;
		bestB1=bigBnow1;

		smallestsic2=sic2+sum(bigBnow2^=0)#bicpen2; 
		sic2=smallestsic2;
		bestB2=bigBnow2;
	end; 
	else do; 
		sic1=sic1+sum(bigBnow1^=0)#bicpen; 
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

/******************************* Use tuning parameter selected above to get graphs *******************************/
bigBnow1=bestB1; bigBnow2=bestB2;
print lambda1 lambda2;


print bigBnow1 bigBnow2; 
	
/* Use the estimated graph structure to get unpenalized estimated regression coefficients */
resids1to2=j(obssize2[<>], pall, 0); resids2to1=j(obssize1[<>], pall, 0); 
resids1to1=j(obssize1[<>], pall, 0); resids2to2=j(obssize2[<>], pall, 0); 
err12=j(obssize1[<>], pall, 0); err21=j(obssize2[<>], pall, 0); 
do j=1 to pall;
	whichrow1=loc(mark1[, j]=1); whichrow2=loc(mark2[, j]=1);
	
	/* first assuming G1 structure */
	whichcol=loc(bigBnow1[, j]^=0); 
	d=ncol(whichcol);
	
	if (d>0) then do;
		/*~~~~~~~~~~~~ center but not standardize w and y ~~~~~~~~~~~~*/
		
		yw=allx2[whichrow2, j]||allx2[whichrow2, whichcol];
		yw=yw-yw[:,]; 
		y2=yw[, 1]; w2=yw[, 2:(d+1)];
		resids1to2[1:obssize2[j], j]=y2-w2*ginv(w2`*w2)*w2`*y2;
		
		yw=allx1[whichrow1, j]||allx1[whichrow1, whichcol];
		yw=yw-yw[:,]; 
		y1=yw[, 1]; w1=yw[, 2:(d+1)];
/*		estB=ginv(w1`*w1)*w1`*y1;*/
/*		resids1to1[1:obssize1[j], j]=y1-w1*estB;*/
		*err21[1:obssize2[j], j]=y2-w2*estB;  /* Use X2 data but B1-hat (obtainedestimating B1 using X1) to find residuals */
		resids1to1[1:obssize1[j], j]=y1-w1*bigBnow1[whichcol, j];
		err21[1:obssize2[j], j]=y2-w2*bigBnow1[whichcol, j];
	end;
	else do; 
		y=allx2[whichrow2, j];
		resids1to2[1:obssize2[j], j]=y-y[:];
		err21[1:obssize2[j], j]=y-y[:]; 

		y=allx1[whichrow1, j];
		resids1to1[1:obssize1[j], j]=y-y[:];
	end;	

	/* second assume G2 structure */
	whichcol=loc(bigBnow2[, j]^=0); 
	d=ncol(whichcol);
	
	if (d>0) then do;
		/*~~~~~~~~~~~~ center but not standardize w and y ~~~~~~~~~~~~*/
		yw=allx1[whichrow1, j]||allx1[whichrow1, whichcol];
		yw=yw-yw[:,]; 
		y1=yw[, 1]; w1=yw[, 2:(d+1)];
		resids2to1[1:obssize1[j], j]=y1-w1*ginv(w1`*w1)*w1`*y1;
		
		yw=allx2[whichrow2, j]||allx2[whichrow2, whichcol];
		yw=yw-yw[:,]; 
		y2=yw[, 1]; w2=yw[, 2:(d+1)];
/*		estB=ginv(w2`*w2)*w2`*y2;*/
/*		resids2to2[1:obssize2[j], j]=y2-w2*estB;*/
/*		err12[1:obssize1[j], j]=y1-w1*estB;*/
		resids2to2[1:obssize2[j], j]=y2-w2*bigBnow2[whichcol, j];
		err12[1:obssize1[j], j]=y1-w1*bigBnow2[whichcol, j];
	end;
	else do; 
		y=allx1[whichrow1, j];
		resids2to1[1:obssize1[j], j]=y-y[:];
		err12[1:obssize1[j], j]=y-y[:]; 

		y=allx2[whichrow2, j];
		resids2to2[1:obssize2[j], j]=y-y[:];
	end;				
end;

/*~~~~~~~~~~~~~~ The following are strategies based on prediction invariance. ~~~~~~~~~~~~~~~~~~*/

/* randomly split resids into two halves corresponding to two sets of experimental conditions under which node j is not intervened */
pvs1to2=j(pall, 3, 1); pvs2to1=j(pall, 3, 1); 
pvs1to1=j(pall, 3, 1); pvs2to2=j(pall, 3, 1); 
discern=j(pall, 1, 1); discern2=j(pall, 1, 1);
PIscore=0; tophowmany=0;
do sp=1 to split;
	do j=1 to pall;
		samsize1=floor(obssize2[j]/2); samsize2=obssize2[j]-samsize1;
		set=j(obssize2[j], 1, 2); 
		cond_idx_inj=cond2_idx[loc(mark2[, j]=1)];
		nocond=ncol(unique(cond_idx_inj));
		if nocond>1 then do;
			firsthalf=floor(nocond/2); 
			call randseed(134+sp+j);
			twoset=ranperm(unique(cond_idx_inj));

			do expcond=1 to firsthalf;
				thiscond=twoset[expcond];
				set[loc(cond_idx_inj=thiscond)]=1; 
			end;
			
			/* Going from G1 to G2 */
			sample1=resids1to2[loc(set=1), j]; sample2=resids1to2[loc(set=2), j]; 
			howcomp=meanvartest(sample1, sample2, samsize1, samsize2);
			pvs1to2[j, ]=howcomp;
			/* Going from G2 to G2 */
		    sample1=resids2to2[loc(set=1), j]; sample2=resids2to2[loc(set=2), j]; 
			howcomp=meanvartest(sample1, sample2, samsize1, samsize2);
			pvs2to2[j, ]=howcomp;
		end;
		else do;
			pvs1to2[j, ]=888; pvs2to2[j, ]=888;
		end;


		samsize1=floor(obssize1[j]/2); samsize2=obssize1[j]-samsize1;
		set=j(obssize1[j], 1, 2); 
		cond_idx_inj=cond1_idx[loc(mark1[, j]=1)];
		nocond=ncol(unique(cond_idx_inj));
		if nocond>1 then do;
			firsthalf=floor(nocond/2); 
			call randseed(13499+sp+j);
			twoset=ranperm(unique(cond_idx_inj));

			do expcond=1 to firsthalf;
				thiscond=twoset[expcond];
				set[loc(cond_idx_inj=thiscond)]=1; 
			end;
			
			/* Going from G2 to G1 */
			sample1=resids2to1[loc(set=1), j]; sample2=resids2to1[loc(set=2), j]; 
			howcomp=meanvartest(sample1, sample2, samsize1, samsize2);
			pvs2to1[j, ]=howcomp;
			/* Going from G1 to G1 */
		    sample1=resids1to1[loc(set=1), j]; sample2=resids1to1[loc(set=2), j]; 
			howcomp=meanvartest(sample1, sample2, samsize1, samsize2);
			pvs1to1[j, ]=howcomp;
		end;
		else do;
			pvs2to1[j, ]=888; pvs1to1[j, ]=888;
		end;

		/* Computer DISCERN score */
		if sp=1 then do;
			commondeno=resids1to1[,j]`*resids1to1[,j]+resids2to2[,j]`*resids2to2[,j];
			discern[j]=(err12[, j]`*err12[, j]+err21[, j]`*err21[, j])/commondeno;
			discern2[j]=(resids2to1[,j]`*resids2to1[,j]+resids1to2[,j]`*resids1to2[,j])/commondeno;
		end; 
	end;
	print sp pvs1to2 pvs1to1 pvs2to1 pvs2to2 discern discern2;
	PIscore=PIscore+(exp(-pvs1to2[, 3])<>exp(-pvs2to1[, 3]))/(exp(-(pvs2to1[, 1:2]<>pvs1to1[, 1:2]))
	<>exp(-(pvs1to2[, 1:2]<>pvs2to2[, 1:2])));
/*	PIscore=PIscore+(exp(-pvs1to2[, 3])<>exp(-pvs2to1[, 3]))/(exp(-pvs1to1[, 3])<>exp(-pvs2to2[, 3]));*/
/*	tophowmany=tophowmany+(sum(pvs1to2[, 3]<pvcutoff))<>(sum(pvs2to1[, 3]<pvcutoff));*/
	tophowmany=tophowmany+sum( ( (pvs2to1[, 1:2]<pvs1to1[, 1:2])#(pvs1to2[, 1:2]<pvs2to2[, 1:2])#
								 (pvs2to1[, 1:2]<pvcutoff)#(pvs1to2[, 1:2]<pvcutoff) )[, #]);
end;
PIscore=PIscore/split;
tophowmany=ceil(tophowmany/split);
print PIscore tophowmany;

if tophowmany>0 then do; 
	keeptrackPI=do(1, pall, 1)`||PIscore;
	call sort(keeptrackPI, 2, 2);
	mydiff=keeptrackPI[1:tophowmany, 1];
	print "PI identifies the differential nodes:" mydiff;

	keeptrackdisc=do(1, pall, 1)`||discern;
	call sort(keeptrackdisc, 2, 2); 
	compdiff=keeptrackdisc[1:tophowmany, 1];
	print "DISERCN-PI identifies the differential nodes:" compdiff;

	keeptrackdisc=do(1, pall, 1)`||discern2;
	call sort(keeptrackdisc, 2, 2); 
	compdiff2=keeptrackdisc[1:tophowmany, 1];
	print "DISERCN2-PI identifies the differential nodes:" compdiff2;

	mydriv_PIpv=FindDriver(mydiff);
	compdriv_PIpv=FindDriver(compdiff);
	comp2driv_PIpv=FindDriver(compdiff2);
	if mydriv_PIpv=-1 then do; 
		print "PI does not identify any driver node.";
	end;
	else do; 
		print "PI identifies driver nodes:" mydriv_PIpv; 
	end;

	if compdriv_PIpv=-1 then do; 
		print "DISCERN-PI does not identify any driver node.";
	end;
	else do; 
		print "DISCERN-PI identifies driver nodes:" compdriv_PIpv; 
	end;

	if comp2driv_PIpv=-1 then do; 
		print "DISCERN2-PI does not identify any driver node.";
	end;
	else do; 
		print "DISCERN2-PI identifies driver nodes:" comp2driv_PIpv; 
	end;
end;
else do; 
	print "PI does not identify any differential node.";
end; 

estpv=j(pall, 2, 0); 
do perm=1 to permno;
	mixdata=allx1[1:intervene_cond1[1], ]//allx2[1:intervene_cond2[1], ];
	call randseed(98765+j+perm);
	unitidx=ranperm(intervene_cond1[1]+intervene_cond2[1]);
	allpermx1=mixdata[unitidx[1:intervene_cond1[1]], ];
	allpermx2=mixdata[unitidx[(intervene_cond1[1]+1):(intervene_cond1[1]+intervene_cond2[1])], ]; 

	do itv=2 to exp_cond;
		mixdata=allx1[(stacknj1[itv-1]+1):stacknj1[itv], ]//allx2[(stacknj2[itv-1]+1):stacknj2[itv], ];
		call randseed(65+j+perm+itv);
		unitidx=ranperm(intervene_cond1[itv]+intervene_cond2[itv]);
		allpermx1=allpermx1//mixdata[unitidx[1:intervene_cond1[itv]], ];
		allpermx2=allpermx2//mixdata[unitidx[(intervene_cond1[itv]+1):(intervene_cond1[itv]+intervene_cond2[itv])], ];
	end;	

	

	bigBnow1_perm=findgraph2(lambda1, allpermx1, bigBini1, mark1); 
	pvmat=getpv2(bigBnow1_perm, allpermx1, mark1);
	bigBnow1_perm=topsort(bigBnow1_perm, pvmat);
	
	bigBnow2_perm=findgraph2(lambda2, allpermx2, bigBini2, mark2); 
	pvmat=getpv2(bigBnow2_perm, allpermx2, mark2);
	bigBnow2_perm=topsort(bigBnow2_perm, pvmat);
	

	do j=1 to pall; 
		whichrow1=loc(mark1[, j]=1); whichrow2=loc(mark2[, j]=1);
			
		permx1=allpermx1[whichrow1, ]; 
		permx2=allpermx2[whichrow2, ];
	

		whichcol=loc(bigBnow1_perm[, j]^=0); 

		d=ncol(whichcol);
		if (d>0) then do;
			/*~~~~~~~~~~~~ center but not standardize w and y ~~~~~~~~~~~~*/
			yw=permx2[, j]||permx2[, whichcol];
			yw=yw-yw[:,]; 
			y2=yw[, 1]; w2=yw[, 2:(d+1)];
			trouble1=ginv(w2`*w2); trouble2=w2`*y2;
			presids1to2=y2-w2*trouble1*trouble2;
			

			yw=permx1[, j]||permx1[, whichcol];
			yw=yw-yw[:,]; 
			y1=yw[, 1]; w1=yw[, 2:(d+1)];
/*			trouble1=ginv(w1`*w1); trouble2=w1`*y1;*/
/*			estB=trouble1*trouble2;*/
/*			presids1to1=y1-w1*estB;*/
/*			perr21=y2-w2*estB; */
			presids1to1=y1-w1*bigBnow1_perm[whichcol, j];
			perr21=y2-w2*bigBnow1_perm[whichcol, j]; 
		end;
		else do; 
			y=permx2[, j];
			presids1to2=y-y[:];
			perr21=y-y[:]; 

			y=permx1[, j];
			presids1to1=y-y[:];
		end;
		
		/* second assume G2 structure */
		whichcol=loc(bigBnow2_perm[, j]^=0); 

		d=ncol(whichcol);
		if (d>0) then do;
			/*~~~~~~~~~~~~ center but not standardize w and y ~~~~~~~~~~~~*/
			yw=permx1[, j]||permx1[, whichcol];
			yw=yw-yw[:,]; 
			y1=yw[, 1]; w1=yw[, 2:(d+1)];
			trouble1=ginv(w1`*w1); trouble2=w1`*y1;
			presids2to1=y1-w1*trouble1*trouble2;
			

			yw=permx2[, j]||permx2[, whichcol];
			yw=yw-yw[:,]; 
			y2=yw[, 1]; w2=yw[, 2:(d+1)];
/*			trouble1=ginv(w2`*w2); trouble2=w2`*y2;*/
/*			estB=trouble1*trouble2;*/
/*			presids2to2=y2-w2*estB;*/
/*			perr12=y1-w1*estB; */
			presids2to2=y2-w2*bigBnow2_perm[whichcol, j];
			perr12=y1-w1*bigBnow1_perm[whichcol, j]; 
		end;
		else do; 
			y=permx1[, j];
			presids2to1=y-y[:];
			perr12=y-y[:]; 

			y=permx2[, j];
			presids2to2=y-y[:];
		end;
		commondeno=presids1to1`*presids1to1+presids2to2`*presids2to2;
		discernj=(perr12`*perr12+perr21`*perr21)/commondeno;
		discern2j=(presids2to1`*presids2to1+presids1to2`*presids1to2)/commondeno;
		estpv[j, 1]=estpv[j, 1]+(discernj>discern[j]);
		estpv[j, 2]=estpv[j, 2]+(discern2j>discern2[j]);
	end;
end;
estpv=estpv/permno; 

pdiscern=BHpv(estpv[, 1]);
if pdiscern=-1 then do;
	print "DISCERN-PM identifies no differential node."; 			
end;
else do;
	print "DISCERN-PM identifies differential nodes:" pdiscern; 
	compdriv_pmpv=FindDriver(pdiscern);
	if compdriv_pmpv=-1 then do;
		print "DISCERN-PM does not identify any driver node";
	end;
	else do;
		print "DISCERN-PM identifies driver nodes:" compdriv_pmpv;
	end;
end;

pdiscern2=BHpv(estpv[, 2]);
if pdiscern2=-1 then do;
	print "DISCERN2-PM identifies no differential node."; 			
end;
else do;
	print "DISCERN2-PM	identifies differential nodes:" pdiscern2; 
	comp2driv_pmpv=FindDriver(pdiscern2);
	if comp2driv_pmpv=-1 then do;
		print "DISCERN2-PM does not identify any driver node";
	end;
	else do;
		print "DISCERN2-PM identifies driver nodes:" comp2driv_pmpv;
	end;
end;

quit;
 
data track;
    t=time(); d=today();
    put 'Time: ' t time9.;
    put 'Date: ' d date9.;
run;

quit;
