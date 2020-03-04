package de.unijena.bioinf.ms.frontend.subtools.decomp;

import de.unijena.bioinf.ms.frontend.subtools.Provide;
import de.unijena.bioinf.ms.frontend.subtools.RootOptions;
import de.unijena.bioinf.ms.frontend.subtools.StandaloneTool;
import de.unijena.bioinf.ms.properties.ParameterConfig;
import picocli.CommandLine;

import java.nio.file.Path;

@CommandLine.Command(name = "decomp", aliases = {}, description = "<STANDALONE> Small tool to decompose masses with given deviation, ionization, chemical alphabet and chemical filter", defaultValueProvider = Provide.Defaults.class, versionProvider = Provide.Versions.class, mixinStandardHelpOptions = true)
public class DecompOptions implements StandaloneTool<DecompWorkflow> {

    @CommandLine.Option(names = {"--ppm", "-p"}, defaultValue = "20", description = "relative mass error in ppm")
    int ppm;

    @CommandLine.Option(names = {"--abs", "-a"}, defaultValue = "0.001", description = "absolute mass error in Dalton")
    double absDeviation;

    @CommandLine.Option(names = {"--nofilter", "-n"}, description = "if set, the molecular formulas are not filtered by their chemical properties")
    boolean noRDBE;

    @CommandLine.Option(names = {"--filter", "-f"}, defaultValue = "NONE", description = "set the strictness of the chemical filter. Allowed levels are STRICT < COMMON < PERMISSIVE < RDBE < NONE.")
    FilterLevel level;

    @CommandLine.Option(names = {"--errors", "-r"}, description = "print mass errors")
    boolean massErrors;

    @CommandLine.Option(names = {"--elements", "-e"}, defaultValue = "CHNOPS", description = "elements which should be used in decomposition. Use syntax 'CH[min-max]N[min-]O[-max]P[num]'")
    private void setAlphabet(String elements) {
        alphabet = new AlphabetParser(elements);
    }
    AlphabetParser alphabet;

    @CommandLine.Option(names = {"--parent"}, description = "If set, the decompositions have to be subformulas of the given parent formula.")
    String parentFormula;

    @CommandLine.Option(names = {"--ion", "-i"}, description = "Ionization mode, e.g. [M+H]+")
    String ionization;

    @CommandLine.Option(names = {"--max-decomps","-d"}, description = "Maximum number of decompositions to be computed. The output might contain less than [--max-decomps] decompositions if there are not more available.")
    Integer maxDecomps;

    @CommandLine.Option(names = {"--mass", "-m"}, description = "Masses that will be decomposed", split = ",")
    double[] masses;

    @CommandLine.Option(names = {"--output", "-o"}, description = "File to with output to.")
    Path out;

    @Override
    public DecompWorkflow makeWorkflow(RootOptions<?,?,?> rootOptions, ParameterConfig config) {
        return new DecompWorkflow(this, rootOptions.getInput());
    }
}
