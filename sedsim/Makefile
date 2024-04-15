#PY = python3
PY = python
SRC = ../../src
REPS = 50

RCFG = 1e5
ROTOR_RADIUS = 10
SPIN_TIME = 900
VISCOSITY = 2

NUMMOLECULES = 30000

all: gen-all

gen-all: gen-mol gen-nuc

test:
#	$(PY) sedsim_test.py test015.test_crosslinking_translation
	$(PY) sedsim_test.py

tln: gen-mol-tln gen-mol-tln-nogene


test-file:
	$(PY) sedsim-driver.py --experiment-file ./sed-control-test.txt --summary-only --summary-out ./results/summary-test.txt

full:
	$(PY) sedsim-driver.py --experiment-file ./sed-control-full.txt --summary-only --summary-out ./results/summary-full-hs.txt

gen-mol-tln:
	$(PY) sedsim-driver.py --template-file ./scer-orf-tpm-occ-length.txt --template-id-field ORF \
	--template-abundance-field Total.TPM.mean --template-length-field LengthTxEst --template-translation-field Occ.mean --number-of-molecules 30000 \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps $(REPS) \
	--experiment-file ./sed-control-mol-tln.txt --summary-only \
	--summary-out ./results/summary-mol-tln-g$(RCFG)-r$(REPS).txt
gen-mol-tln-nogene:
	$(PY) sedsim-driver.py --template-file ./scer-orf-tpm-occ-length.txt --template-id-field ORF \
	--template-abundance-field Total.TPM.mean --template-length-field LengthTxEst --number-of-molecules 30000 \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps $(REPS) \
	--experiment-file ./sed-control-mol-tln.txt --summary-only \
	--summary-out ./results/summary-mol-tln-nogene-g$(RCFG)-r$(REPS).txt
gen-tln:
	$(PY) sedsim-driver.py --template-file ./scer-orf-tpm-occ-length.txt --template-id-field ORF \
	--template-abundance-field Total.TPM.mean --template-length-field LengthTxEst --template-translation-field Occ.mean --number-of-molecules 30000 \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps $(REPS) \
	--experiment-file ./sed-control-tln.txt \
	--out ./results/sedctrl-tln-g$(RCFG)-r$(REPS).txt --summary-out ./results/summary-tln-g$(RCFG)-r$(REPS).txt
gen-tln-nogene:
	$(PY) sedsim-driver.py --template-file ./scer-orf-tpm-occ-length.txt --template-id-field ORF \
	--template-abundance-field Total.TPM.mean --template-length-field LengthTxEst --number-of-molecules 30000 \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps $(REPS) \
	--experiment-file ./sed-control-tln.txt \
	--out ./results/sedctrl-tln-nogene-g$(RCFG)-r$(REPS).txt --summary-out ./results/summary-tln-nogene-g$(RCFG)-r$(REPS).txt
gen-tln-test:
	$(PY) sedsim-driver.py --template-file ./scer-orf-tpm-occ-length.txt --template-id-field ORF \
	--template-abundance-field Total.TPM.mean --template-length-field LengthTxEst --template-translation-field Occ.mean --number-of-molecules 300000 \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps 1 \
	--experiment-file ./sed-control-tln.txt --summary-only \
	--out ./results/sedctrl-tln-test-g$(RCFG)-r$(REPS).txt --summary-out ./results/summary-tln-test-g$(RCFG)-r$(REPS).txt
gen-tln-nogene-test-summary:
	$(PY) sedsim-driver.py --template-file ./scer-orf-tpm-occ-length.txt --template-id-field ORF \
	--template-abundance-field Total.TPM.mean --template-length-field LengthTxEst --number-of-molecules 30000 \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps 1 \
	--experiment-file ./sed-control-tln.txt --summary-only \
	--summary-out ./results/summary-tln-nogene-test-g$(RCFG)-r1.txt


gen-nuc:
	$(PY) sedsim-driver.py --template-file $(SRC)/scer-orf-mrna-length.txt --template-id-field ORF \
	--template-abundance-field mrna --template-length-field LengthTxEst --number-of-molecules $(NUMMOLECULES) \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps $(REPS) \
	--experiment-file ./sed-control-nuc.txt \
	--out ./results/sedctrl-nuc-g$(RCFG)-r$(REPS).txt --summary-out ./results/summary-nuc-g$(RCFG)-r$(REPS).txt
gen-mol:
	$(PY) sedsim-driver.py --template-file $(SRC)/scer-orf-mrna-length.txt --template-id-field ORF \
	--template-abundance-field mrna --template-length-field LengthTxEst --number-of-molecules $(NUMMOLECULES) \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps $(REPS) \
	--experiment-file ./sed-control-mol.txt \
	--out ./results/sedctrl-mol-g$(RCFG)-r$(REPS).txt --summary-out ./results/summary-mol-g$(RCFG)-r$(REPS).txt
gen-dom:
	$(PY) sedsim-driver.py --template-file $(SRC)/scer-orf-mrna-length.txt --template-id-field ORF \
	--template-abundance-field mrna --template-length-field LengthTxEst --number-of-molecules $(NUMMOLECULES) \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps $(REPS) \
	--experiment-file ./sed-control-mol-dominant.txt \
	--out ./results/sedctrl-dom-g$(RCFG)-r$(REPS).txt --summary-out ./results/summary-dom-g$(RCFG)-r$(REPS).txt

gen-full: gen-nuc-full gen-mol-full 
#gen-dom-full
FULL_REPS = 25
gen-nuc-full-old:
	$(PY) sedsim-driver.py --template-file $(SRC)/scer-orf-mrna-length.txt --template-id-field ORF \
	--template-abundance-field mrna --template-length-field LengthTxEst --number-of-molecules $(NUMMOLECULES) \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps $(FULL_REPS) \
	--experiment-file ./sed-control-nuc.txt --summary-only \
	--summary-out ./results/summary-nuc-g$(RCFG)-full-r$(FULL_REPS).txt

gen-mol-full:
	$(PY) sedsim-driver.py --experiment-file ./sed-control-mol.txt --summary-only --summary-out ./results/summary-mol-full.txt
gen-nuc-full:
	$(PY) sedsim-driver.py --experiment-file ./sed-control-nuc.txt --summary-only --summary-out ./results/summary-nuc-full.txt
gen-mol-full-old:
	$(PY) sedsim-driver.py --template-file $(SRC)/scer-orf-mrna-length.txt --template-id-field ORF \
	--template-abundance-field mrna --template-length-field LengthTxEst --number-of-molecules $(NUMMOLECULES) \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps $(FULL_REPS) \
	--experiment-file ./sed-control-mol.txt --summary-only \
	--summary-out ./results/summary-mol-g$(RCFG)-full-r$(FULL_REPS).txt
gen-dom-full:
	$(PY) sedsim-driver.py --template-file $(SRC)/scer-orf-mrna-length.txt --template-id-field ORF \
	--template-abundance-field mrna --template-length-field LengthTxEst --number-of-molecules $(NUMMOLECULES) \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps $(FULL_REPS) \
	--experiment-file ./sed-control-mol-dominant.txt --summary-only \
	--summary-out ./results/summary-dom-g$(RCFG)-full-r$(FULL_REPS).txt

# Testing.
sum-mol:
	$(PY) sedsim-driver.py --template-file $(SRC)/scer-orf-mrna-length.txt --template-id-field ORF \
	--template-abundance-field mrna --template-length-field LengthTxEst --number-of-molecules $(NUMMOLECULES) \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps 2 \
	--experiment-file ./sed-control-mol.txt \
	--summary-only --summary-out ./results/summary-mol-g$(RCFG)-r2.txt

gen-test:
	$(PY) sedsim-driver.py --template-file ./scer-orf-tpm-occ-length.txt --template-id-field ORF \
	--template-abundance-field Total.TPM.mean --template-length-field LengthTxEst --template-translation-field Occ.mean --number-of-molecules 30000 \
	--rcf-g $(RCFG) --spin-time $(SPIN_TIME) --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps 2 \
	--experiment-file ./sed-control-test.txt \
	--out ./results/test.txt --summary-out ./results/test-summary.txt

calibration:
	$(PY) sedsim-driver.py --template-file ./calibration-sizes.txt --template-id-field id \
	--template-abundance-field abund --template-length-field length --number-of-molecules 10000 \
	--rcf-g $(RCFG) --spin-time 900 --rotor-radius $(ROTOR_RADIUS) --viscosity $(VISCOSITY) --reps 100 \
	--experiment-file ./sed-control-calibration.txt \
	--out ./results/calibration-g$(RCFG).txt --summary-out ./results/calibration-summary-g$(RCFG).txt

