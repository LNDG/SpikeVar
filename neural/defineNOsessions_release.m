%
% List of all NO sessions that are part of the public data release.
% [please note that many sessions don't exist here because they belong to different experiments]
% 
% urut&acarlson&mfaraut/Jul17
function [NOsessions, NO_listOf_allUsable] = defineNOsessions_release()

NOsessions=[];  % struct of all sessions

%P9S1 same day
c=5; 
NOsessions(c).session='P9HMH_032306'; 
NOsessions(c).sessionID='P9S1';
NOsessions(c).EXPERIMENTIDLearn=80; % Learning block trials in eventRaw.mat file
NOsessions(c).EXPERIMENTIDRecog=81; % Memory test block trials in eventRaw.mat file
NOsessions(c).taskDescr='NO';   % folder name in which data are stored (in 'events' and 'sorted' folders)
NOsessions(c).variant=1;            % 3 variants of the experiment corresponding to 3 sets of stimuli (cf. Table 4)
NOsessions(c).blockIDLearn=1;       % used to load the stimuli numbers for the Learning block
NOsessions(c).blockIDRecog=2;       % used to load the stimuli numbers for the Memory block
NOsessions(c).patient_nb=1;         % patient number
NOsessions(c).patientsession=1;     % rank of the session for a given patient (if a same patient did the task more than ones)
NOsessions(c).diagnosisCode=1;     % location of epileptic focal point: 0)not localized 1)Right Mesial Temporal 2)Left Mesial Temporal 3)Right Neocortical Temporal 4)Left Neocortical Temporal 5)Right Lateral Frontal 6)Left Lateral Frontal 7)Bilateral Independent Temporal 8)Bilateral Independent Frontal 9)Right Other 10)Left Other

%P9S3. same day.
c=6;
NOsessions(c).session='P9HMH_032506';
NOsessions(c).sessionID='P9S3';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=0; 
NOsessions(c).patient_nb=1;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=1;

%P10S2 same day
c=7;
NOsessions(c).session='P10HMH_092206';
NOsessions(c).sessionID='P10S2';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=0; 
NOsessions(c).patient_nb=2;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=6;

%P11S1 NO same day
c=9;
NOsessions(c).session='P11HMH_110906';
NOsessions(c).sessionID='P11S1';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=3;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=5;

%P14S1. same day.
c=17;
NOsessions(c).session='P14HMH_062107';
NOsessions(c).sessionID='P14S1';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=1; 
NOsessions(c).patient_nb=4;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=7;

%P14S2. same day. (reduced)
c=18;
NOsessions(c).session='P14HMH_062307';
NOsessions(c).sessionID='P14S2';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=1; 
NOsessions(c).patient_nb=4;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=7;

%P15S1. same day.
c=20;
NOsessions(c).session='P15HMH_091307';
NOsessions(c).sessionID='P15S1';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=0; 
NOsessions(c).patient_nb=5;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

%P15S2. same day.
c=21;
NOsessions(c).session='P15HMH_091407';
NOsessions(c).sessionID='P15S2';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=5;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=1;

%P16S1. same day.
c=23;
NOsessions(c).session='P16HMH_101207';
NOsessions(c).sessionID='P16S1';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=1; 
NOsessions(c).patient_nb=6;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=5;

%P16S2. same day.
c=24;
NOsessions(c).session='P16HMH_101507';
NOsessions(c).sessionID='P16S2';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=1; %sort on raw, LFP on rawLFP. 
NOsessions(c).patient_nb=6;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=5;

%P18S2. same day. (reduced)
c=26;
NOsessions(c).session='P18HMH_061608';
NOsessions(c).sessionID='P18S2';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=1; 
NOsessions(c).patient_nb=7;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

%P19S1. same day. (reduced)
c=27;
NOsessions(c).session='P19HMH_062608';
NOsessions(c).sessionID='P19S1';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=1; 
NOsessions(c).patient_nb=8;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=6;

%P19S2. same day. 
c=28;
NOsessions(c).session='P19HMH_062708';
NOsessions(c).sessionID='P19S2';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=1; 
NOsessions(c).patient_nb=8;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=6;

%P17S1. same day. (reduced)
c=32;
NOsessions(c).session='P17HMH_052208';
NOsessions(c).sessionID='P17S1';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).phaseShiftCorrection=1; 
NOsessions(c).patient_nb=9;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=6;

%P21S1 NO2 same day
c=38;
NOsessions(c).session='P21HMH_012209';
NOsessions(c).sessionID='P21S1';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;
NOsessions(c).taskDescr='NO2';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=10;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=0;

%P21S1 NO1 same day
c=39;
NOsessions(c).session='P21HMH_012209';
NOsessions(c).sessionID='P21S1';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO1';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=10;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=0;

%P21S2 NO3 same day 
c=41;
NOsessions(c).session='P21HMH_012309';
NOsessions(c).sessionID='P21S2';
NOsessions(c).EXPERIMENTIDLearn=88;
NOsessions(c).EXPERIMENTIDRecog=89;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=10;
NOsessions(c).patientsession=3;
NOsessions(c).diagnosisCode=0;

%P23S2 NO3 same day
c=43;
NOsessions(c).session='P23HMH_022109';
NOsessions(c).sessionID='P23S2';
NOsessions(c).EXPERIMENTIDLearn=88;
NOsessions(c).EXPERIMENTIDRecog=89;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=11;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=2;

%P23S2 NO2 same day (reduced -- 50 learning trials)
c=44;
NOsessions(c).session='P23HMH_022109';
NOsessions(c).sessionID='P23S2';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=11;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=2;

%P23S4 NO1 same day
c=47;
NOsessions(c).session='P23HMH_022309';
NOsessions(c).sessionID='P23S4';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=11;
NOsessions(c).patientsession=3;
NOsessions(c).diagnosisCode=2;

%P28 HMH, 031310
c=48;
NOsessions(c).session='P28HMH';
NOsessions(c).sessionID='P28S2';
NOsessions(c).EXPERIMENTIDLearn=88;
NOsessions(c).EXPERIMENTIDRecog=89;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=12;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

%P29 HMH, 121810
c=49;
NOsessions(c).session='P29HMH';
NOsessions(c).sessionID='P29S2';
NOsessions(c).EXPERIMENTIDLearn=88;
NOsessions(c).EXPERIMENTIDRecog=89;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=13;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=10;

%P31 HMH, 042611
c=50;
NOsessions(c).session='P31HMH';
NOsessions(c).sessionID='P31';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=14;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

%P33 HMH, 072011
c=52;
NOsessions(c).session='P33HMH';
NOsessions(c).sessionID='P33S2';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=15;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=2;

%P42 HMH, 033013
c=54;
NOsessions(c).session='P42HMH';
NOsessions(c).sessionID='P42';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  %Also has 82, same day
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=16;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=0;

%P42 HMH, 033013; second block (same day)
c=55;
NOsessions(c).session='P42HMH';
NOsessions(c).sessionID='P42';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=82;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=3; 
NOsessions(c).patient_nb=16;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=0;

%P43S1 HMH, 041213
c=56;
NOsessions(c).session='P43HMH_S2';
NOsessions(c).sessionID='P43S1';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=17;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=2;

%P27 (KF), 022310
c=58;
NOsessions(c).session='P27HMH_022310';
NOsessions(c).sessionID='P27S2';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=18;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=7;

%P24CS, session1
c=59;
NOsessions(c).session='P24CS_091912';
NOsessions(c).sessionID='P24CSS2';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=19;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=0;

%P24CS, session2
c=60;
NOsessions(c).session='P24CS_092112';
NOsessions(c).sessionID='P24CSS2';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=19;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=0;

%P25CS, session1
c=61;
NOsessions(c).session='P25CS_092712';
NOsessions(c).sessionID='P25CSS1';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=20;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=7;

%P25CS, session2
c=63;
NOsessions(c).session='P25CS_092812';
NOsessions(c).sessionID='P25CSS2';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=20;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=7;

%P26CS, session1
c=64;
NOsessions(c).session='P26CS_121112';
NOsessions(c).sessionID='P26CSS1';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=21;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=2;

%P26CS, session3
c=66;
NOsessions(c).session='P26CS_121412';
NOsessions(c).sessionID='P26CSS3';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=21;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=2;

%P27CS, session1
c=67;
NOsessions(c).session='P27CS_011913';
NOsessions(c).sessionID='P27CSS1';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=22;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=2;

%P44 HMH, 
c=68;
NOsessions(c).session='P44HMH_s3';
NOsessions(c).sessionID='P44HMHs3';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=23;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

%P29CS, s1
c=69;
NOsessions(c).session='P29CS_103013';
NOsessions(c).sessionID='P29CSs1';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=24;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=4;

%P29CS, s2, with eye tracking
c=70;
NOsessions(c).session='P29CS_103113';
NOsessions(c).sessionID='P29CSs2';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=24;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=4;

%var3, no eye track
c=73;
NOsessions(c).session='P31CS_020314';
NOsessions(c).sessionID='P31CSs2';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=25;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=4;

%var2, no eye track
c=74;
NOsessions(c).session='P32CS_021314';
NOsessions(c).sessionID='P32CSs1';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=26;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=0;

%var3, no eye track
c=76;
NOsessions(c).session='P33CS_032714';
NOsessions(c).sessionID='P33CS';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;  
NOsessions(c).taskDescr='NO3';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=27;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

%var2, no eye track
c=77;
NOsessions(c).session='P33CS_032714';
NOsessions(c).sessionID='P33CS';
NOsessions(c).EXPERIMENTIDLearn=81;
NOsessions(c).EXPERIMENTIDRecog=82;  
NOsessions(c).taskDescr='NO2';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=27;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=1;

%var1, same day, with ET
c=78;
NOsessions(c).session='P33CS_033014';
NOsessions(c).sessionID='P33CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=27;
NOsessions(c).patientsession=3;
NOsessions(c).diagnosisCode=1;

%NOvar3, with ET
c=85;
NOsessions(c).session='P34CS_121414';
NOsessions(c).sessionID='P34CS';
NOsessions(c).EXPERIMENTIDLearn=83;
NOsessions(c).EXPERIMENTIDRecog=84;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=28;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=7;

%var1, NO P47HMH, no ET
c=92;
NOsessions(c).session='P47HMH_062014';
NOsessions(c).sessionID='P47HMH';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=29;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

%NO P39CS, var1, with eye track
c=93;
NOsessions(c).session='P39CS_081515';
NOsessions(c).sessionID='P39CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=30;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=9;

%NO P37CS, var3, no eye track
c=96;
NOsessions(c).session='P37CS_032515';
NOsessions(c).sessionID='P37CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=31;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

% NO P47HMH, var3 no ET
c=97;
NOsessions(c).session='P47HMH_062214';
NOsessions(c).sessionID='P47HMH';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=29;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=1;


% NO P48 HMH var 2
c = 98;
NOsessions(c).session='P48HMH_122014';
NOsessions(c).sessionID='P48HMH';
NOsessions(c).EXPERIMENTIDLearn=82;
NOsessions(c).EXPERIMENTIDRecog=83;  
NOsessions(c).taskDescr='NO2';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=32;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=2;

% NO P48 HMH var 1
c=99;
NOsessions(c).session='P48HMH_122014';
NOsessions(c).sessionID='P48HMH';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO1';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=32;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=2;


%NO P39CS, var3
c=100;
NOsessions(c).session='P39CS_080515';
NOsessions(c).sessionID='P39CS';
NOsessions(c).EXPERIMENTIDLearn=82;
NOsessions(c).EXPERIMENTIDRecog=83;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=30;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=9;

%NO P40CS, var3
c=101;
NOsessions(c).session='P40CS_010916';
NOsessions(c).sessionID='P40CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=33;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=9;

%NO P38CS, var3
c=102;
NOsessions(c).session='P38CS_070115';
NOsessions(c).sessionID='P38CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=34;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

%NO P51HMH, var3
c=104;
NOsessions(c).session='P51HMH_021516';
NOsessions(c).sessionID='P51HMH';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=35;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=8;

%NO P51HMH, var 1
c=105;
NOsessions(c).session='P51HMH_021416';
NOsessions(c).sessionID='P51HMH';
NOsessions(c).EXPERIMENTIDLearn=82;
NOsessions(c).EXPERIMENTIDRecog=83;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2; 
NOsessions(c).patient_nb=35;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=8;

%NO P42CS, var1
c=111;
NOsessions(c).session='P42CS_081416';
NOsessions(c).sessionID='P42CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=36;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=0;

%NO P42CS, var2
c=112;
NOsessions(c).session='P42CS_081516';
NOsessions(c).sessionID='P42CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=2;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=36;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=0;

%NO P43CS, var3
c=113;
NOsessions(c).session='P43CS_110316';
NOsessions(c).sessionID='P43CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=37;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=2;

%NO P44CS, var3
c=114;
NOsessions(c).session='P44CS_090516';
NOsessions(c).sessionID='P44CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=38;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

%NO P47CS, var3
c=115;
NOsessions(c).session='P47CS_022017';
NOsessions(c).sessionID='P47CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=39;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=1;

%NO P48CS, var3
c=116;
NOsessions(c).session='P48CS_031017';
NOsessions(c).sessionID='P48CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=40;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=2;

%NO P49CS, var3
c=117;
NOsessions(c).session='P49CS_052217';
NOsessions(c).sessionID='P49CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=41;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=2;

%NO P49CS, var1
c=118;
NOsessions(c).session='P49CS_052617';
NOsessions(c).sessionID='P49CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=41;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=2;

%NO P51CS, var3
c=119;
NOsessions(c).session='P51CS_070117';
NOsessions(c).sessionID='P51CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO3';
NOsessions(c).variant=3;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=42;
NOsessions(c).patientsession=1;
NOsessions(c).diagnosisCode=0;

%NO P51CS, var1
c=120;
NOsessions(c).session='P51CS_070117';
NOsessions(c).sessionID='P51CS';
NOsessions(c).EXPERIMENTIDLearn=80;
NOsessions(c).EXPERIMENTIDRecog=81;  
NOsessions(c).taskDescr='NO1';
NOsessions(c).variant=1;
NOsessions(c).blockIDLearn=1; 
NOsessions(c).blockIDRecog=2;
NOsessions(c).patient_nb=42;
NOsessions(c).patientsession=2;
NOsessions(c).diagnosisCode=0;

%=====
%
NO_listOf_allUsable = [];
usable=[];
for k=1:length(NOsessions)
    if ~isempty( NOsessions(k).session )
        NO_listOf_allUsable = [NO_listOf_allUsable k];
    end
end
