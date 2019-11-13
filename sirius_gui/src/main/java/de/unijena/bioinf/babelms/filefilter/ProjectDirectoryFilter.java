package de.unijena.bioinf.babelms.filefilter;

import javax.swing.filechooser.FileFilter;
import java.io.File;

public class ProjectDirectoryFilter extends FileFilter {
    @Override
    public boolean accept(File f) {
       return f.isDirectory();
    }

    @Override
    public String getDescription() {
        return "Sirius project directory.";
    }
}