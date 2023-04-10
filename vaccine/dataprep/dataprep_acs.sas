/* Prep ACS data	*/
/* 3/11/2023 RA		*/

options nodate nocenter nonumber ls=256 formdlim=" " formchar="|----|+|---+=|-/\<>*";

libname acs "./data_ACS2019";

proc format cntlin=acs.usa_00002_f library=work;
run;
data acs(keep=YEAR CLUSTER STRATA PERWT REPWTP1-REPWTP80
			statefips state male educcat educcat_coll educ_: racecat race_: hisp racethcat raceth_: agecat age_: fb_:
			incomecat income_: incomemisscat incomemiss_: income6cat income6_: mscat mscat_coll ms_:);
set acs.Usa_00002;
*state fips code;
statefips = STATEFIP;
state = fipstate(statefips);
*binary gender; *no missing;
/*1 Male*/
/*2 Female*/
male = (SEX = 1);
*education;
if EDUCD in(2,11,12,14,15,16,17,22,23,25,26) then educcat = 1; *< HS; 
	else if EDUCD in(30,40,50,61) then educcat = 2; *some HS;
	else if EDUCD in(63,64) then educcat = 3; *HS/GED;
	else if EDUCD in(65,71) then educcat = 4; *some coll no degree;
	else if EDUCD in(81) then educcat = 5; *assoc;
	else if EDUCD in(101) then educcat = 6; *bachelors;
	else if EDUCD in(114,115,116) then educcat = 7; *grad/prof;
educ_lths = (educcat = 1);
educ_somehs = (educcat = 2);
educ_hs = (educcat = 3);
educ_somecollnodeg = (educcat = 4);
educ_assoc = (educcat = 5);
educ_bach = (educcat = 6);
educ_graddeg = (educcat = 7);
*education levels collapsed;
educ_lthsgrad = (educcat in(1,2));
educ_hsless = (educcat in(1,2,3));
educ_somecoll = (educcat in(4,5));
educcat_coll = 1*educ_hsless + 2*educ_somecoll + 3*educ_bach + 4*educ_graddeg;
*education for Facebook;
fb_educcat_coll = 1*educ_lthsgrad + 2*educ_hs + 3*educ_somecoll + 4*educ_bach + 5*educ_graddeg;
*race;
/*1 White*/
/*2 Black/African American*/
/*3 American Indian or Alaska Native*/
/*4 Chinese*/
/*5 Japanese*/
/*6 Other Asian or Pacific Islander*/
/*7 Other race, nec*/
/*8 Two major races*/
/*9 Three or more major races*/
if RACE = 1 then racecat = 1; *white alone;
	else if RACE = 2 then racecat = 2; *black alone;
	else if RACE in(4,5,6) then racecat = 3; *asian alone (incl pacific islander);
	else if RACE in(3,7,8,9) then racecat = 4; *other incl AIAN/more than one race;
race_white = (racecat = 1);
race_black = (racecat = 2);
race_asian = (racecat = 3);
race_other = (racecat = 4);
*ethnicity;
/*0 Not Hispanic*/
/*1 Mexican*/
/*2 Puerto Rican*/
/*3 Cuban*/
/*4 Other*/
/*9 Not Reported*/
hisp = (HISPAN in(1,2,3,4));
*race/ethnicity combined;
if racecat = 1 and hisp = 0 then racethcat = 1; *NH white;
	else if racecat = 2 and hisp = 0 then racethcat = 2; *NH Black;
	else if RACE in(3,4,5,6,7,8,9) and hisp = 0 then racethcat = 3; *NH other;
	else if hisp = 1 then racethcat = 4; *Hispanic;
if racethcat > . then do;
	raceth_whitenh = (racethcat = 1);
	raceth_blacknh = (racethcat = 2);
	raceth_othernh = (racethcat = 3);
	raceth_hisp = (racethcat = 4);
end;
*race/ethnicity combined for Facebook;
if RACE = 1 and hisp = 0 then fb_racethcat = 1; *NH white;
	else if RACE = 2 and hisp = 0 then fb_racethcat = 2; *NH Black;
	else if RACE = 3 and hisp = 0 then fb_racethcat = 3; *NH AIAN;
	else if RACE in(4,5,6) and hisp = 0 then fb_racethcat = 4; *NH Asian (incl pacific islander);
	else if RACE in(7,8,9) and hisp = 0 then fb_racethcat = 5; *NH Other incl multiple races;
	else if hisp = 1 then fb_racethcat = 6; *Hispanic;
if fb_racethcat > . then do;
	fb_raceth_whitenh = (fb_racethcat = 1);
	fb_raceth_blacknh = (fb_racethcat = 2);
	fb_raceth_aiannh = (fb_racethcat = 3);
	fb_raceth_asianaapinh = (fb_racethcat = 4);
	fb_raceth_othernh = (fb_racethcat = 5);
	fb_raceth_hisp = (fb_racethcat = 6);
end;
*age;
age_18_29 = (18 <= AGE <= 29);
age_30_39 = (30 <= AGE <= 39);
age_40_49 = (40 <= AGE <= 49);
age_50_59 = (50 <= AGE <= 59);
age_60_69 = (60 <= AGE <= 69);
age_70plus = (AGE >= 70);
agecat = 1*age_18_29 + 2*age_30_39 + 3*age_40_49 + 4*age_50_59 + 5*age_60_69 + 6*age_70plus;
*age for Facebook (different categories);
fb_age_18_24 = (18 <= AGE <= 24);
fb_age_25_34 = (25 <= AGE <= 34);
fb_age_35_44 = (35 <= AGE <= 44);
fb_age_45_54 = (45 <= AGE <= 54);
fb_age_55_64 = (55 <= AGE <= 64);
fb_age_65_74 = (65 <= AGE <= 74);
fb_age_75plus = (AGE >= 75);
fb_agecat = 1*fb_age_18_24 + 2*fb_age_25_34 + 3*fb_age_35_44 + 4*fb_age_45_54 + 5*fb_age_55_64 + 6*fb_age_65_74 + 7*fb_age_75plus;
*income;
if FTOTINC not in(9999999) then do;
	income_1 = (FTOTINC <= 24999);
	income_2 = (25000 <= FTOTINC <= 34999);
	income_3 = (35000 <= FTOTINC <= 49999);
	income_4 = (50000 <= FTOTINC <= 74999);
	income_5 = (75000 <= FTOTINC <= 99999);
	income_6 = (100000 <= FTOTINC <= 149999);
	income_7 = (150000 <= FTOTINC <= 199999);
	income_8 = (FTOTINC >= 200000);
end;
incomecat = 1*income_1 + 2*income_2 + 3*income_3 + 4*income_4 + 5*income_5 + 6*income_6 + 7*income_7 + 8*income_8; 
*include missing income as a cateogory;
incomemisscat = incomecat;
if incomemisscat = . then incomemisscat = 0; *missing income;
incomemiss_0 = (FTOTINC = 9999999);
incomemiss_1 = (FTOTINC <= 24999);
incomemiss_2 = (25000 <= FTOTINC <= 34999);
incomemiss_3 = (35000 <= FTOTINC <= 49999);
incomemiss_4 = (50000 <= FTOTINC <= 74999);
incomemiss_5 = (75000 <= FTOTINC <= 99999);
incomemiss_6 = (100000 <= FTOTINC <= 149999);
incomemiss_7 = (150000 <= FTOTINC <= 199999);
incomemiss_8 = (FTOTINC >= 200000);
*6-level income;
if FTOTINC not in(9999999) then do;
	income6_1 = (FTOTINC <= 24999);
	income6_2 = (25000 <= FTOTINC <= 49999);
	income6_3 = (50000 <= FTOTINC <= 74999);
	income6_4 = (75000 <= FTOTINC <= 99999);
	income6_5 = (100000 <= FTOTINC <= 149999);
	income6_6 = (FTOTINC >= 150000);
end;
income6cat = 1*income6_1 + 2*income6_2 + 3*income6_3 + 4*income6_4 + 5*income6_5 + 6*income6_6;
*marital status; *no missing data;
ms_married = (MARST in(1,2));
ms_widowed = (MARST = 5);
ms_divorced = (MARST = 4);
ms_separated = (MARST = 3);
ms_single = (MARST = 6);
mscat = 1*ms_married + 2*ms_widowed + 3*ms_divorced + 4*ms_separated + 5*ms_single;
*marital status collapsed;
ms_divsep = (MARST in(3,4));
mscat_coll = 1*ms_married + 2*ms_widowed + 3*ms_divsep + 4*ms_single;
run;

/* Covariance matrices for PPMA */

/**** for HPS - WITH INCOME ****/
	/* Get weighted SSCP matrix */
	proc reg data=acs outsscp=sscp;
	var male
		educ_somecoll educ_bach educ_graddeg /* ref=educ_hslss */
		race_black race_asian race_other /* ref=race_white */
		hisp
		age_18_29 age_30_39 age_40_49 age_50_59 age_60_69 /* ref=age_70plus */
		income_1 income_2 income_3 income_4 income_5 income_6 income_7 /* ref=income_8 */
	;
	weight PERWT;
	run;quit;
	data sscp;
	set sscp(where=(_TYPE_="SSCP"));
	run;
	/* Turn into covariance matrix */
	proc princomp data=sscp(TYPE=SSCP) outstat=covmat cov;
	var male
		educ_somecoll educ_bach educ_graddeg /* ref=educ_hslss */
		race_black race_asian race_other /* ref=race_white */
		hisp
		age_18_29 age_30_39 age_40_49 age_50_59 age_60_69 /* ref=age_70plus */
		income_1 income_2 income_3 income_4 income_5 income_6 income_7 /* ref=income_8 */
	;
	run;
	data covmat;
	set covmat;
	if _TYPE_ in ("SCORE", "EIGENVAL","N") then delete;
	run;
	proc export data=covmat outfile="./data_ACS2019/acs_covmat_hps_with_income.csv" DBMS=CSV replace;
	run;
	proc datasets nolist;
	delete sscp covmat;
	quit;

  /**** for HPS - NO INCOME ****/
	/* Get weighted SSCP matrix */
	proc reg data=acs outsscp=sscp;
	var male
		educ_somecoll educ_bach educ_graddeg /* ref=educ_hslss */
		race_black race_asian race_other /* ref=race_white */
		hisp
		age_18_29 age_30_39 age_40_49 age_50_59 age_60_69 /* ref=age_70plus */
	;
	weight PERWT;
	run;quit;
	data sscp;
	set sscp(where=(_TYPE_="SSCP"));
	run;
	/* Turn into covariance matrix */
	proc princomp data=sscp(TYPE=SSCP) outstat=covmat cov;
	var male
		educ_somecoll educ_bach educ_graddeg /* ref=educ_hslss */
		race_black race_asian race_other /* ref=race_white */
		hisp
		age_18_29 age_30_39 age_40_49 age_50_59 age_60_69 /* ref=age_70plus */
	;
	run;
	data covmat;
	set covmat;
	if _TYPE_ in ("SCORE", "EIGENVAL","N") then delete;
	run;
	proc export data=covmat outfile="./data_ACS2019/acs_covmat_hps.csv" DBMS=CSV replace;
	run;
	proc datasets nolist;
	delete sscp covmat;
	quit;

/**** for AXIOS-IPSOS ****/
	/* Get weighted SSCP matrix */
	proc reg data=acs outsscp=sscp;
	var male
		educ_somecoll educ_bach educ_graddeg /* ref=educ_hsless */
		raceth_blacknh raceth_othernh raceth_hisp /* ref=raceth_whitenh */
		age_18_29 age_30_39 age_40_49 age_50_59 age_60_69 /* ref=age_70plus */
		income6_1 income6_2 income6_3 income6_4 income6_5 /* ref=income6_6 */
	;
	weight PERWT;
	run;quit;
	data sscp;
	set sscp(where=(_TYPE_="SSCP"));
	run;
	/* Turn into covariance matrix */
	proc princomp data=sscp(TYPE=SSCP) outstat=covmat cov;
	var male
		educ_somecoll educ_bach educ_graddeg /* ref=educ_hsless */
		raceth_blacknh raceth_othernh raceth_hisp /* ref=raceth_whitenh */
		age_18_29 age_30_39 age_40_49 age_50_59 age_60_69 /* ref=age_70plus */
		income6_1 income6_2 income6_3 income6_4 income6_5 /* ref=income6_6 */
	;
	run;
	data covmat;
	set covmat;
	if _TYPE_ in ("SCORE", "EIGENVAL","N") then delete;
	run;
	proc export data=covmat outfile="./data_ACS2019/acs_covmat_axios.csv" DBMS=CSV replace;
	run;
/**** for FACEBOOK ****/
	/* Get weighted SSCP matrix */
	proc reg data=acs outsscp=sscp;
		var male
		educ_hs educ_somecoll educ_bach educ_graddeg /* ref=educ_lthsgrad */
		fb_raceth_blacknh fb_raceth_aiannh fb_raceth_asianaapinh fb_raceth_othernh fb_raceth_hisp /* ref=fb_raceth_whitenh */
		fb_age_18_24 fb_age_25_34 fb_age_35_44 fb_age_45_54 fb_age_55_64 fb_age_65_74 /* ref=fb_age_75plus */
	;
	weight PERWT;
	run;quit;
	data sscp;
	set sscp(where=(_TYPE_="SSCP"));
	run;
	/* Turn into covariance matrix */
	proc princomp data=sscp(TYPE=SSCP) outstat=covmat cov;
		var male
		educ_hs educ_somecoll educ_bach educ_graddeg /* ref=educ_lthsgrad */
		fb_raceth_blacknh fb_raceth_aiannh fb_raceth_asianaapinh fb_raceth_othernh fb_raceth_hisp /* ref=fb_raceth_whitenh */
		fb_age_18_24 fb_age_25_34 fb_age_35_44 fb_age_45_54 fb_age_55_64 fb_age_65_74 /* ref=fb_age_75plus */
	;
	run;
	data covmat;
	set covmat;
	if _TYPE_ in ("SCORE", "EIGENVAL","N") then delete;
	run;
	proc export data=covmat outfile="./data_ACS2019/acs_covmat_facebook.csv" DBMS=CSV replace;
	run;

/**** for FACEBOOK - AGE*GENDER ONLY ****/
	data acs;
	set acs;
	male_fb_age_18_24 = male*fb_age_18_24;
	male_fb_age_25_34 = male*fb_age_25_34;
	male_fb_age_35_44 = male*fb_age_35_44;
	male_fb_age_45_54 = male*fb_age_45_54;
	male_fb_age_55_64 = male*fb_age_55_64;
	male_fb_age_65_74 = male*fb_age_65_74;
	run;
	/* Get weighted SSCP matrix */
	proc reg data=acs outsscp=sscp;
	var male
		fb_age_18_24 fb_age_25_34 fb_age_35_44 fb_age_45_54 fb_age_55_64 fb_age_65_74 /* ref=fb_age_75plus */
		male_fb_age_18_24 male_fb_age_25_34 male_fb_age_35_44 male_fb_age_45_54 male_fb_age_55_64 male_fb_age_65_74
	;
	weight PERWT;
	run;quit;
	data sscp;
	set sscp(where=(_TYPE_="SSCP"));
	run;
	/* Turn into covariance matrix */
	proc princomp data=sscp(TYPE=SSCP) outstat=covmat cov;
		var male
		fb_age_18_24 fb_age_25_34 fb_age_35_44 fb_age_45_54 fb_age_55_64 fb_age_65_74 /* ref=fb_age_75plus */
		male_fb_age_18_24 male_fb_age_25_34 male_fb_age_35_44 male_fb_age_45_54 male_fb_age_55_64 male_fb_age_65_74
	;
	run;
	data covmat;
	set covmat;
	if _TYPE_ in ("SCORE", "EIGENVAL","N") then delete;
	run;
	proc export data=covmat outfile="./data_ACS2019/acs_covmat_facebook_AGE_GENDER.csv" DBMS=CSV replace;
	run;

