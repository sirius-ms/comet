package de.unijena.bioinf.MassDecomposer.cli;


import de.unijena.bioinf.ChemistryBase.chem.Element;
import de.unijena.bioinf.ChemistryBase.chem.PeriodicTable;
import de.unijena.bioinf.MassDecomposer.Alphabet;
import de.unijena.bioinf.MassDecomposer.Chemistry.ChemicalAlphabet;
import de.unijena.bioinf.MassDecomposer.Interval;
import de.unijena.bioinf.MassDecomposer.ValencyAlphabet;

import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class AlphabetParser {

    private final HashMap<Element, Interval> boundaryMap;
    private final ChemicalAlphabet alphabet;


    private final static Pattern REGEXP = Pattern.compile("([A-Z][a-z]*)(\\[(\\d*)(-\\d*)?\\])?");

    /*
     * Syntax: C[1-10]H[-3]N[1-]O[2-6]P[1-9]S[5]
     */
    public AlphabetParser(String s) {
        final PeriodicTable table = PeriodicTable.getInstance();
        final Matcher matcher = REGEXP.matcher(s);
        final HashMap<Element, Interval> map = new HashMap<Element, Interval>();
        while (matcher.find()) {
            final String elementName = matcher.group(1);
            final Element element = table.getByName(elementName);
            if (element == null) throw new RuntimeException("Unknown element with name '" + elementName + "'");
            if (map.containsKey(element)) {
                throw new RuntimeException("Element '" + elementName + "' is contained twice in given alphabet.");
            }
            final Interval boundary;
            if (matcher.group(2) != null && !matcher.group(2).isEmpty()) {
                final String begin = matcher.group(3);
                final String end = matcher.group(4);
                if (begin != null && !begin.isEmpty()) {
                    final int a = Integer.parseInt(begin);
                    if (end != null && end.length() > 1) {
                        boundary = new Interval(a, Integer.parseInt(end.substring(1)));
                    } else if (end.equals("-")) {
                        boundary = new Interval(a, Integer.MAX_VALUE);
                    } else {
                        boundary = new Interval(a, a);
                    }
                } else if (end != null && end.length() > 1) {
                    boundary = new Interval(0, Integer.parseInt(end.substring(1)));
                } else {
                    throw new RuntimeException("Cannot parse boundary '" + elementName + matcher.group(2)  +"'");
                }
            } else {
                boundary = new Interval(0, Integer.MAX_VALUE);
            }
            map.put(element, boundary);
        }
        this.alphabet = new ChemicalAlphabet(map.keySet().toArray(new Element[map.size()]));
        this.boundaryMap = map;
    }

    public HashMap<Element, Interval> getBoundary() {
        return boundaryMap;
    }

    public ChemicalAlphabet getAlphabet() {
        return alphabet;
    }

}
