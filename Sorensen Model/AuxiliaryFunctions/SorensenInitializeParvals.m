%========================================================================================================================================
% SCRIPT SorensenInitializeParvals.m: assignment of initial values to named parameters
% Sorensen V01.01.41 20190724 (Gemini 13.01.06, BMLib 10.0.2, Autocoder 02.11.09, coded 24-Jul-2019 16:06:54)
%========================================================================================================================================


global PARMIN PARMAX PARDETM;

nPars = 142;

PARMIN = zeros(142,1);
PARMAX = zeros(142,1);
PARDETM = zeros(142,1);
PARMIN(001) = -30; PARMAX(001) = 1440; PARDETM(001)=0;
PARMIN(002) = 0; PARMAX(002) = 10080; PARDETM(002)=0;
PARMIN(003) = 0; PARMAX(003) = 600; PARDETM(003)=0;
PARMIN(004) = 0; PARMAX(004) = 10; PARDETM(004)=0;
PARMIN(005) = 0; PARMAX(005) = 10; PARDETM(005)=0;
PARMIN(006) = 0; PARMAX(006) = 10; PARDETM(006)=0;
PARMIN(007) = 0.01; PARMAX(007) = 50; PARDETM(007)=0;
PARMIN(008) = 0.1; PARMAX(008) = 20; PARDETM(008)=0;
PARMIN(009) = 0; PARMAX(009) = 10; PARDETM(009)=0;
PARMIN(010) = 0; PARMAX(010) = 10; PARDETM(010)=0;
PARMIN(011) = 0; PARMAX(011) = 10; PARDETM(011)=0;
PARMIN(012) = 0; PARMAX(012) = 10; PARDETM(012)=0;
PARMIN(013) = 0; PARMAX(013) = 10; PARDETM(013)=0;
PARMIN(014) = 0; PARMAX(014) = 10; PARDETM(014)=0;
PARMIN(015) = 0; PARMAX(015) = 10; PARDETM(015)=0;
PARMIN(016) = 0; PARMAX(016) = 10; PARDETM(016)=0;
PARMIN(017) = 0; PARMAX(017) = 10; PARDETM(017)=0;
PARMIN(018) = 0; PARMAX(018) = 10; PARDETM(018)=0;
PARMIN(019) = 0; PARMAX(019) = 10; PARDETM(019)=0;
PARMIN(020) = 0; PARMAX(020) = 10; PARDETM(020)=0;
PARMIN(021) = 0; PARMAX(021) = 10; PARDETM(021)=0;
PARMIN(022) = 0; PARMAX(022) = 10; PARDETM(022)=0;
PARMIN(023) = 0; PARMAX(023) = 10; PARDETM(023)=0;
PARMIN(024) = 0.01; PARMAX(024) = 50; PARDETM(024)=0;
PARMIN(025) = 0; PARMAX(025) = 10; PARDETM(025)=0;
PARMIN(026) = 0; PARMAX(026) = 100; PARDETM(026)=0;
PARMIN(027) = 0; PARMAX(027) = 100; PARDETM(027)=0;
PARMIN(028) = 0; PARMAX(028) = 100; PARDETM(028)=0;
PARMIN(029) = 0; PARMAX(029) = 100; PARDETM(029)=0;
PARMIN(030) = 0; PARMAX(030) = 100; PARDETM(030)=0;
PARMIN(031) = 0; PARMAX(031) = 100; PARDETM(031)=0;
PARMIN(032) = 0; PARMAX(032) = 100; PARDETM(032)=0;
PARMIN(033) = 0; PARMAX(033) = 100; PARDETM(033)=0;
PARMIN(034) = 0; PARMAX(034) = 100; PARDETM(034)=0;
PARMIN(035) = 0; PARMAX(035) = 100; PARDETM(035)=0;
PARMIN(036) = 0; PARMAX(036) = 100; PARDETM(036)=0;
PARMIN(037) = 0; PARMAX(037) = 100; PARDETM(037)=0;
PARMIN(038) = 0; PARMAX(038) = 100; PARDETM(038)=0;
PARMIN(039) = 0; PARMAX(039) = 100; PARDETM(039)=0;
PARMIN(040) = 0; PARMAX(040) = 100; PARDETM(040)=0;
PARMIN(041) = 0; PARMAX(041) = 100; PARDETM(041)=0;
PARMIN(042) = 0; PARMAX(042) = 20; PARDETM(042)=0;
PARMIN(043) = 0; PARMAX(043) = 100; PARDETM(043)=0;
PARMIN(044) = 0; PARMAX(044) = 100; PARDETM(044)=0;
PARMIN(045) = 0; PARMAX(045) = 100; PARDETM(045)=0;
PARMIN(046) = 0; PARMAX(046) = 100; PARDETM(046)=0;
PARMIN(047) = 0; PARMAX(047) = 100; PARDETM(047)=0;
PARMIN(048) = 0; PARMAX(048) = 100; PARDETM(048)=0;
PARMIN(049) = 0; PARMAX(049) = 20; PARDETM(049)=0;
PARMIN(050) = 0; PARMAX(050) = 100; PARDETM(050)=0;
PARMIN(051) = 0; PARMAX(051) = 100; PARDETM(051)=0;
PARMIN(052) = 0; PARMAX(052) = 100; PARDETM(052)=0;
PARMIN(053) = 0; PARMAX(053) = 100; PARDETM(053)=0;
PARMIN(054) = 0; PARMAX(054) = 100; PARDETM(054)=0;
PARMIN(055) = 0; PARMAX(055) = 100; PARDETM(055)=0;
PARMIN(056) = 0; PARMAX(056) = 10; PARDETM(056)=0;
PARMIN(057) = 0; PARMAX(057) = 10; PARDETM(057)=0;
PARMIN(058) = 0; PARMAX(058) = 10; PARDETM(058)=0;
PARMIN(059) = 0; PARMAX(059) = 10; PARDETM(059)=0;
PARMIN(060) = 0; PARMAX(060) = 10; PARDETM(060)=0;
PARMIN(061) = 0; PARMAX(061) = 10; PARDETM(061)=0;
PARMIN(062) = 0; PARMAX(062) = 10; PARDETM(062)=0;
PARMIN(063) = 0; PARMAX(063) = 10; PARDETM(063)=0;
PARMIN(064) = 0; PARMAX(064) = 10; PARDETM(064)=0;
PARMIN(065) = 0; PARMAX(065) = 10; PARDETM(065)=0;
PARMIN(066) = 0; PARMAX(066) = 10; PARDETM(066)=0;
PARMIN(067) = 0; PARMAX(067) = 10; PARDETM(067)=0;
PARMIN(068) = 0; PARMAX(068) = 10; PARDETM(068)=0;
PARMIN(069) = 0; PARMAX(069) = 10; PARDETM(069)=0;
PARMIN(070) = 0; PARMAX(070) = 10; PARDETM(070)=0;
PARMIN(071) = 0; PARMAX(071) = 50; PARDETM(071)=0;
PARMIN(072) = 0; PARMAX(072) = 10; PARDETM(072)=0;
PARMIN(073) = 0; PARMAX(073) = 10; PARDETM(073)=0;
PARMIN(074) = 0; PARMAX(074) = 10; PARDETM(074)=0;
PARMIN(075) = 0; PARMAX(075) = 10; PARDETM(075)=0;
PARMIN(076) = 0; PARMAX(076) = 10; PARDETM(076)=0;
PARMIN(077) = 0; PARMAX(077) = 10; PARDETM(077)=0;
PARMIN(078) = 0; PARMAX(078) = 1; PARDETM(078)=0;
PARMIN(079) = 0; PARMAX(079) = 100000; PARDETM(079)=0;
PARMIN(080) = 0; PARMAX(080) = 1e+06; PARDETM(080)=0;
PARMIN(081) = 0; PARMAX(081) = 1; PARDETM(081)=0;
PARMIN(082) = 0; PARMAX(082) = 10; PARDETM(082)=0;
PARMIN(083) = 0; PARMAX(083) = 1; PARDETM(083)=0;
PARMIN(084) = 0; PARMAX(084) = 1; PARDETM(084)=0;
PARMIN(085) = 0.1; PARMAX(085) = 1000; PARDETM(085)=0;
PARMIN(086) = 0.1; PARMAX(086) = 1000; PARDETM(086)=0;
PARMIN(087) = 0; PARMAX(087) = 10; PARDETM(087)=0;
PARMIN(088) = 0.1; PARMAX(088) = 1000; PARDETM(088)=0;
PARMIN(089) = 0; PARMAX(089) = 100; PARDETM(089)=0;
PARMIN(090) = 0; PARMAX(090) = 100; PARDETM(090)=0;
PARMIN(091) = 0; PARMAX(091) = 100; PARDETM(091)=0;
PARMIN(092) = 0; PARMAX(092) = 100; PARDETM(092)=0;
PARMIN(093) = 0; PARMAX(093) = 100; PARDETM(093)=0;
PARMIN(094) = 0; PARMAX(094) = 100; PARDETM(094)=0;
PARMIN(095) = 0; PARMAX(095) = 100; PARDETM(095)=0;
PARMIN(096) = 0; PARMAX(096) = 100; PARDETM(096)=0;
PARMIN(097) = 0; PARMAX(097) = 100; PARDETM(097)=0;
PARMIN(098) = 0; PARMAX(098) = 150000; PARDETM(098)=0;
PARMIN(099) = 0; PARMAX(099) = 150000; PARDETM(099)=0;
PARMIN(100) = -10; PARMAX(100) = 100; PARDETM(100)=0;
PARMIN(101) = -10; PARMAX(101) = 100; PARDETM(101)=0;
PARMIN(102) = 0; PARMAX(102) = 20; PARDETM(102)=0;
PARMIN(103) = 0; PARMAX(103) = 80000; PARDETM(103)=0;
PARMIN(104) = -3; PARMAX(104) = 1440; PARDETM(104)=0;
PARMIN(105) = -3; PARMAX(105) = 1440; PARDETM(105)=0;
PARMIN(106) = 0; PARMAX(106) = 1000; PARDETM(106)=1;
PARMIN(107) = 0; PARMAX(107) = 1000; PARDETM(107)=1;
PARMIN(108) = 0; PARMAX(108) = 1000; PARDETM(108)=1;
PARMIN(109) = 0; PARMAX(109) = 1000; PARDETM(109)=1;
PARMIN(110) = 0; PARMAX(110) = 1000; PARDETM(110)=1;
PARMIN(111) = 0; PARMAX(111) = 1000; PARDETM(111)=1;
PARMIN(112) = 0; PARMAX(112) = 1000; PARDETM(112)=1;
PARMIN(113) = 0; PARMAX(113) = 1000; PARDETM(113)=1;
PARMIN(114) = 0; PARMAX(114) = 1; PARDETM(114)=1;
PARMIN(115) = 0; PARMAX(115) = 1; PARDETM(115)=1;
PARMIN(116) = 0; PARMAX(116) = 1; PARDETM(116)=1;
PARMIN(117) = 0; PARMAX(117) = 1; PARDETM(117)=1;
PARMIN(118) = 0; PARMAX(118) = 200100; PARDETM(118)=1;
PARMIN(119) = 0; PARMAX(119) = 10000; PARDETM(119)=1;
PARMIN(120) = 0.1; PARMAX(120) = 50; PARDETM(120)=1;
PARMIN(121) = 0.1; PARMAX(121) = 50; PARDETM(121)=1;
PARMIN(122) = 0.1; PARMAX(122) = 50; PARDETM(122)=1;
PARMIN(123) = 0.1; PARMAX(123) = 50; PARDETM(123)=1;
PARMIN(124) = 0.1; PARMAX(124) = 50; PARDETM(124)=1;
PARMIN(125) = 0.1; PARMAX(125) = 50; PARDETM(125)=1;
PARMIN(126) = 0.1; PARMAX(126) = 50; PARDETM(126)=1;
PARMIN(127) = 0; PARMAX(127) = 20; PARDETM(127)=1;
PARMIN(128) = 0; PARMAX(128) = 10; PARDETM(128)=1;
PARMIN(129) = 0; PARMAX(129) = 10; PARDETM(129)=1;
PARMIN(130) = 0; PARMAX(130) = 10; PARDETM(130)=1;
PARMIN(131) = 0; PARMAX(131) = 10; PARDETM(131)=1;
PARMIN(132) = 0; PARMAX(132) = 100; PARDETM(132)=1;
PARMIN(133) = 0; PARMAX(133) = 10; PARDETM(133)=1;
PARMIN(134) = 0; PARMAX(134) = 10; PARDETM(134)=1;
PARMIN(135) = 0; PARMAX(135) = 20; PARDETM(135)=1;
PARMIN(136) = 0; PARMAX(136) = 20; PARDETM(136)=1;
PARMIN(137) = 0; PARMAX(137) = 1000; PARDETM(137)=1;
PARMIN(138) = 0; PARMAX(138) = 1000; PARDETM(138)=1;
PARMIN(139) = 0; PARMAX(139) = 100; PARDETM(139)=1;
PARMIN(140) = 0; PARMAX(140) = 10; PARDETM(140)=1;
PARMIN(141) = 0; PARMAX(141) = 100; PARDETM(141)=1;
PARMIN(142) = 0; PARMAX(142) = 100; PARDETM(142)=1;

Tzero                = -30; % (bigtheta 001);
Tend                 = 200; % (bigtheta 002);
Tdelta               = 0.1; % (bigtheta 003);
QfloGB               = 0.59; % (bigtheta 004);
VolGBV               = 0.35; % (bigtheta 005);
VolBI                = 0.45; % (bigtheta 006);
TdifB                = 2.1; % (bigtheta 007);
GlucH0               = 5.07333; % (bigtheta 008);
GammaBGU             = 0.388889; % (bigtheta 009);
QfloGL               = 1.26; % (bigtheta 010);
QfloGK               = 1.01; % (bigtheta 011);
QfloGP               = 1.51; % (bigtheta 012);
QfloGH               = 4.37; % (bigtheta 013);
GammaRBCU            = 0.0555556; % (bigtheta 014);
VolGH                = 1.38; % (bigtheta 015);
QfloGJ               = 1.01; % (bigtheta 016);
VolGJ                = 1.12; % (bigtheta 017);
GammaJGU             = 0.111111; % (bigtheta 018);
QfloGA               = 0.25; % (bigtheta 019);
VolGL                = 2.51; % (bigtheta 020);
VolGK                = 0.66; % (bigtheta 021);
VolGPV               = 1.04; % (bigtheta 022);
VolPI                = 6.74; % (bigtheta 023);
TdifGP               = 5; % (bigtheta 024);
GammaBPGU            = 0.194444; % (bigtheta 025);
beta0PGU             = 7.03; % (bigtheta 026);
beta1PGU             = 6.52; % (bigtheta 027);
beta2PGU             = 0.338; % (bigtheta 028);
beta3PGU             = 5.82; % (bigtheta 029);
beta0HGP             = 2.7; % (bigtheta 030);
beta1HGP             = 0.388852; % (bigtheta 031);
tauCgon              = 65; % (bigtheta 032);
beta2HGP             = 1.21; % (bigtheta 033);
beta3HGP             = 1.14; % (bigtheta 034);
beta4HGP             = 1.66; % (bigtheta 035);
beta5HGP             = 0.887748; % (bigtheta 036);
tauInsu              = 25; % (bigtheta 037);
beta6HGP             = 1.42; % (bigtheta 038);
beta7HGP             = 1.41; % (bigtheta 039);
beta8HGP             = 0.62; % (bigtheta 040);
beta9HGP             = 0.504543; % (bigtheta 041);
GammaHGP0            = 0.861111; % (bigtheta 042);
beta0HGU             = 2; % (bigtheta 043);
beta1HGU             = 0.549306; % (bigtheta 044);
beta2HGU             = 5.66; % (bigtheta 045);
beta3HGU             = 5.66; % (bigtheta 046);
beta4HGU             = 2.44; % (bigtheta 047);
beta5HGU             = 1.4783; % (bigtheta 048);
GammaHGU0            = 0.111111; % (bigtheta 049);
beta0KGE             = 0.394444; % (bigtheta 050);
beta1KGE             = 0.394444; % (bigtheta 051);
beta2KGE             = 0.198; % (bigtheta 052);
beta3KGE             = 25.5556; % (bigtheta 053);
beta4KGE             = 1.834; % (bigtheta 054);
beta5KGE             = 0.0872; % (bigtheta 055);
QfloIB               = 0.45; % (bigtheta 056);
VolIB                = 0.26; % (bigtheta 057);
VolIH                = 0.99; % (bigtheta 058);
QfloIL               = 0.9; % (bigtheta 059);
QfloIK               = 0.72; % (bigtheta 060);
QfloIP               = 1.05; % (bigtheta 061);
QfloIH               = 3.12; % (bigtheta 062);
VolIJ                = 0.94; % (bigtheta 063);
QfloIJ               = 0.72; % (bigtheta 064);
VolIL                = 1.14; % (bigtheta 065);
QfloIA               = 0.18; % (bigtheta 066);
FracLIC              = 0.4; % (bigtheta 067);
FracKIC              = 0.3; % (bigtheta 068);
VolIK                = 0.51; % (bigtheta 069);
VolIPV               = 0.74; % (bigtheta 070);
TdifIP               = 20; % (bigtheta 071);
FracPIC              = 0.15; % (bigtheta 072);
beta1PIR             = 3.27; % (bigtheta 073);
beta2PIR             = 7.33333; % (bigtheta 074);
beta3PIR             = 2.879; % (bigtheta 075);
beta4PIR             = 3.02; % (bigtheta 076);
beta5PIR             = 1.11; % (bigtheta 077);
KappaRinsu           = 0.00794; % (bigtheta 078);
Rinsu0               = 44310; % (bigtheta 079);
KappaRinsuPotn       = 4025; % (bigtheta 080);
KappaPotnPtgt        = 0.0482; % (bigtheta 081);
KappaPinhPrp         = 0.931; % (bigtheta 082);
EMME1                = 0.00747; % (bigtheta 083);
EMME2                = 0.0958; % (bigtheta 084);
InsuPV0              = 91; % (bigtheta 085);
Cgon0                = 11.48; % (bigtheta 086);
GammaMCC             = 0.91; % (bigtheta 087);
VolC                 = 11.31; % (bigtheta 088);
beta0PCR             = 2.93; % (bigtheta 089);
beta1PCR             = 2.1; % (bigtheta 090);
beta2PCR             = 4.18; % (bigtheta 091);
beta3PCR             = 0.621325; % (bigtheta 092);
beta4PCR             = 1.31; % (bigtheta 093);
beta5PCR             = 0.61; % (bigtheta 094);
beta6PCR             = 1.06; % (bigtheta 095);
beta7PCR             = 0.471419; % (bigtheta 096);
Func20               = 0; % (bigtheta 097);
GammaIVG0            = 0; % (bigtheta 098);
GammaIVGin           = 64.81; % (bigtheta 099);
TimeIVG              = -3; % (bigtheta 100);
TimeIVGend           = 0; % (bigtheta 101);
GammaIVI0            = 0; % (bigtheta 102);
GammaIVIin           = 0; % (bigtheta 103);
TimeIVI              = 0; % (bigtheta 104);
TimeIVIend           = 0; % (bigtheta 105);
InsuH0               = 107.059; % (bigtheta 106);
InsuK0               = 74.9412; % (bigtheta 107);
InsuB0               = 107.059; % (bigtheta 108);
InsuJ0               = 107.059; % (bigtheta 109);
InsuPI0              = 40.9651; % (bigtheta 110);
InsuL0               = 151.488; % (bigtheta 111);
GammaBPIR            = 130.879; % (bigtheta 112);
GammaPIC0            = 16.8618; % (bigtheta 113);
Pprp0                = 0.19032; % (bigtheta 114);
Ptgt0                = 0.158572; % (bigtheta 115);
Pinh0                = 0.19032; % (bigtheta 116);
Potn0                = 0.158572; % (bigtheta 117);
InitialRinsu0        = 108507; % (bigtheta 118);
Secr0                = 128.53; % (bigtheta 119);
GlucPV0              = 4.94456; % (bigtheta 120);
GlucK0               = 5.07333; % (bigtheta 121);
GlucBV0              = 4.4142; % (bigtheta 122);
GlucJ0               = 4.96332; % (bigtheta 123);
GlucL0               = 5.58039; % (bigtheta 124);
GlucBI0              = 2.59938; % (bigtheta 125);
GlucPI0              = 4.80032; % (bigtheta 126);
MIPGU0               = 0.992859; % (bigtheta 127);
MCHGP0               = 1; % (bigtheta 128);
MC0HGP0              = 1; % (bigtheta 129);
MIHGP0               = 1; % (bigtheta 130);
MIHGPinf0            = 1; % (bigtheta 131);
MGHGP0               = 1; % (bigtheta 132);
MIHGU0               = 1; % (bigtheta 133);
MIHGUinf0            = 1; % (bigtheta 134);
MGHGU0               = 1; % (bigtheta 135);
GammaKGE0            = 0.000236777; % (bigtheta 136);
GammaLIC0            = 90.8929; % (bigtheta 137);
GammaKIC0            = 23.1247; % (bigtheta 138);
MGPCR0               = 1; % (bigtheta 139);
MIPCR0               = 1; % (bigtheta 140);
GammaPCC0            = 10.4468; % (bigtheta 141);
GammaBPCR            = 10.4468; % (bigtheta 142);

% run DetermineParameters on the named parameters before assigning all values to bigtheta
SorensenDetermineParameters;

% build bigtheta
bigtheta = zeros(nPars,1);
SorensenParvals2Bigtheta;

% check that all parameters are Kosher
errcode = SorensenCheckParameters(bigtheta);
if(errcode > 0)
   error(sprintf('\nERROR specifying parameter %03d: out of bounds',errcode));
end

%========================================================================================================================================