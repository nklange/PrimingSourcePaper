# Repository for "Explaining the association between repetition priming and source memory: No evidence for a contribution of recognition or fluency" (Nicholas Lange, Christopher J. Berry)
## Raw data

Inside _Data_ folder, sorted by experiment.

Exp 1: CIDRS14
Exp 2: CIDRS19
Exp 3: CIDRS21
Exp 4: CIDRS22

For each experiment, participants' raw data is in individual folders.

In CIDRS14, CIDRS19, and CIDRS21:
* StudStim.csv (studied words)
* XX_ID_CHECK.csv (identification phase: target, response, correct)
* dataSTUDY_XX.csv (study phase)
* dataTEST_XX.csv (test phase)
* and header files for these files ([...]Heads.csv)

In CIDRS22:

- StudStim.csv (studied words)
- XX_ID_CHECK.csv (identification phase: target, response, correct)
- dataSTUDY_XX.csv (study phase)
- dataTESTID_XX.csv (test phase identification task)
- dataTESTJUD_XX.csv (test phase source memory task)
- and header files for these files ([...]Heads.csv)

In pre-processing of the raw data, information from these files is combined. Identification of items was manually corrected for typographical errors, i.e., 'duck13' or `ducck' is counted as 'correct' for 'duck'.

## Pre-processing of data

Raw data are pre-processed in R script (StandardAnalysis/CIDS_preprocessing.R) to produce pre-processed files CIDRSXX_fulldatascreened.R

In all columns "-99" indicates a missing value. In CIDRS19 and CIDSR21 this is used to indicate that no recognition task was completed. In CIDRS22, the recognition task columns are omitted.

ExclReason column indicates which trials are not analysed further (with a reason given). The only trials that are analysed are those marked ``valid''.

Columns in _fulldatascreened.R files:
* TestTrial: number of trial at test
* Block: experimental block in test phase
* StimType: Stimulus Type
	* 1 = Old / Presented at study
	* 2 = New
* StimID: ID number of stimulus word
* CIDCorrect: Identification of word in CID task
	* 1 = Correct
	* 2 = Incorrect
* RecJud: Recognition judgment
	* 1 = ``Old''
	* 2 = ``New''
* recRT: Recognition judgment RT during memory rating phase from presentation of the scale to response
* SDTclass: Type of response in signal detection terms
	* 1 = Hit (``Old''|Old)
	* 2 = False Alarm (``Old''|New)
	* 3 = Miss (``New''|Old)
	* 4 = Correct rejection (``New''|New)
* RecConf: Recognition confidence rating 
	* 1 = sure new
	* 2 = probably new
	* 3 = guess new
	* 4 = guess old
	* 5 = probably old
	* 6 = sure old
* SourceJud: Source judgment
	* 1 = ``Top''
	* 2 = ``Bottom''
* SourceRT: Source judgment RT during memory rating phase from presentation of the scale to response
* SourceCorr: Accuracy of source judgment
	* 1 = ``Correct''
	* 2 = ``Incorrect'', 
	* 0 = [response made to new item]
* SourceConf: Source Confidence rating 
	* 1 = sure bottom
	* 2 = probably bottom
	* 3 = guess bottom
	* 4 = guess top
	* 5 = probably top
	* 6 = sure top
* Bottom1Top2New0: Source manipulation
	* 1 = presented at bottom at study
	* 2 = presented at top at study
	* 0 = not presented at study
* SubjID: ID of participant
* SubjAge: Age of participant
* M1F2: Gender of participant
	* 1 = Male
	* 2 = Female
* TotExpDur: Total duration of the experiment
* ExclReason (Reasons for trials to be excluded in pre-processing): 
	* valid = trial is valid
	* NoResponse = no response in CID task
	* CIDincorrect = word was incorrectly identified in CID task
	* TooSlow = Identification RT of CID task exceeded mean + 3SD for this participant
	* TooFast = Identification RT < 200 ms
* item: Item shown at test in CID task
* item_presponse: Participant’s typed word / identification in CID task

## Analysis

Inside *StandardAnalysis* folder, subsections of analysis commented in the .R files

summaryfunctions.R - some general functions
CIDS_preprocessing.R - pre-processing of raw data
CIDRS14_standardanalysis.R - analysis CIDRS14 (Exp 1)
CIDRS19_standardanalysis.R - analysis CIDRS19 (Exp 2)
CIDRS21_standardanalysis.R - analysis CIDRS21 (Exp 3)
CIDRS22_standardanalysis.R - analysis CIDRS22 (Exp 4)

