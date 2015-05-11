/*
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2015 Kai Dührkop
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with SIRIUS.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.unijena.bioinf.FragmentationTreeConstruction.computation.inputValidator;

import de.unijena.bioinf.ChemistryBase.chem.ChemicalAlphabet;
import de.unijena.bioinf.ChemistryBase.chem.ElectronIonization;
import de.unijena.bioinf.ChemistryBase.chem.FormulaConstraints;
import de.unijena.bioinf.ChemistryBase.chem.Ionization;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.ChemistryBase.ms.EIIntensityDeviation;
import de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment;
import de.unijena.bioinf.ChemistryBase.ms.MutableMeasurementProfile;
import de.unijena.bioinf.FragmentationTreeConstruction.model.Ms2ExperimentImpl;

public class GCMSMissingValueValidator extends MissingValueValidator {

    @Override
    public Ms2Experiment validate(Ms2Experiment originalInput, Warning warn, boolean repair) throws InvalidException {
        final Ms2ExperimentImpl input = new Ms2ExperimentImpl(originalInput);
        if (input.getMs1Spectra() == null || input.getMs1Spectra().isEmpty()) throw new InvalidException("Missing MS1 spectra");
        checkIonization(warn, repair, input);
        checkMergedMs1(warn, repair, input);
        checkNeutralMass(warn, repair, input);
        checkMeasurementProfile(warn, repair, input);
        return input;
    }

    @Override
    protected void checkIonization(Warning warn, boolean repair, Ms2ExperimentImpl input) {
        if (input.getIonization() == null) {
            throwOrWarn(warn, repair, "No ionization is given");
            final Ionization ion = new ElectronIonization();
            input.setIonization(ion);
        }
    }

    @Override
    protected void checkNeutralMass(Warning warn, boolean repair, Ms2ExperimentImpl input) {
        if (input.getMoleculeNeutralMass() == 0 || !validDouble(input.getMoleculeNeutralMass(), false)) {
            if (input.getMolecularFormula() != null) {
                throwOrWarn(warn, repair, "Neutral mass is missing, but formula known.");
                input.setMoleculeNeutralMass(input.getMolecularFormula().getMass());
            } else if (input.getIonMass() != 0){
                throwOrWarn(warn, repair, "Neutral mass is missing, but ion mass known.");
                input.setMoleculeNeutralMass(input.getIonization().subtractFromMass(input.getIonMass()));
            }
        }
    }

    @Override
    protected void checkMeasurementProfile(Warning warn, boolean repair, Ms2ExperimentImpl input) {
        if (input.getMeasurementProfile() == null) throw new InvalidException("Measurement profile is missing");
        final MutableMeasurementProfile profile = new MutableMeasurementProfile(input.getMeasurementProfile());
        if (profile.getFormulaConstraints() == null) {
            throwOrWarn(warn, repair, "Measurement profile: Formula constraints are missing");
            profile.setFormulaConstraints(new FormulaConstraints(new ChemicalAlphabet()));
        }
        // get at least one deviation
        Deviation dev = profile.getAllowedMassDeviation();
        if (dev == null) dev = profile.getStandardMs1MassDeviation();
        if (profile.getAllowedMassDeviation() == null) {
            System.out.println("set allowed "+profile.getStandardMs1MassDeviation());

            throwOrWarn(warn, repair && dev!=null, "Measurement profile: Maximal allowed Mass deviation is missing");
            profile.setAllowedMassDeviation(dev);
        }
        if (profile.getStandardMs1MassDeviation() == null) {
            throwOrWarn(warn, repair && dev!=null, "Measurement profile: MS1 deviation is missing");
            profile.setStandardMs1MassDeviation(dev);
        }
        if (profile.getStandardMs2MassDeviation() != null) {
            throwOrWarn(warn, repair && dev!=null, "Measurement profile: standardMs2MassDeviation not null. But EI uses no MS2. StandardMs2MassDeviation is set to same as MS1.");
        }
        profile.setStandardMs2MassDeviation(profile.getStandardMs1MassDeviation());
        if (!dev.getClass().equals(EIIntensityDeviation.class)) throwOrWarn(warn, false, "Deviation is not from class 'EIIntensityDeviation' (but intensity is necessary for calculations).");
        input.setMeasurementProfile(profile);
    }

    private boolean validDouble(double val, boolean mayNegative) {
        return !Double.isInfinite(val) && !Double.isNaN(val) && (mayNegative || val >= 0d);
    }

}
