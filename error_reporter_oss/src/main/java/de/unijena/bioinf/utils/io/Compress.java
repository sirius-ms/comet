/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schilller University.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with SIRIUS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.utils.io;
/**
 * Created by Markus Fleischauer (markus.fleischauer@gmail.com)
 * as part of the sirius
 * 29.09.16.
 */

import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**
 * @author Markus Fleischauer (markus.fleischauer@gmail.com)
 */
public class Compress {
    public static void compressToZipArchive(File zipFile, File... srcFiles) {
        try (
                FileOutputStream fos = new FileOutputStream(zipFile);
                ZipOutputStream zos = new ZipOutputStream(fos)
        ) {
            // create byte buffer
            byte[] buffer = new byte[1024];

            for (File srcFile : srcFiles) {
                try (FileInputStream fis = new FileInputStream(srcFile)) {

                    // begin writing a new ZIP entry, positions the stream to the start of the entry data
                    zos.putNextEntry(new ZipEntry(srcFile.getName()));
                    int length;
                    while ((length = fis.read(buffer)) > 0) {
                        zos.write(buffer, 0, length);
                    }
                    zos.closeEntry();
                    // close the InputStream
                    fis.close();
                } catch (IOException ioe) {
                    LoggerFactory.getLogger(Compress.class).error("Could not Compress " + srcFile.getAbsolutePath(), ioe);
                    throw ioe;
                }
            }
            // close the ZipOutputStream
            zos.close();
        } catch (IOException ioe) {
            LoggerFactory.getLogger(Compress.class).error("Could not Create zip archive " + zipFile.getAbsolutePath(), ioe);
        }
    }

    public static void compressToZipArchive(File zipFile, Map<InputStream, String> input) {
        try {
            compressToZipArchive(new FileOutputStream(zipFile), input);
        } catch (FileNotFoundException e) {
            LoggerFactory.getLogger(Compress.class).error("Could not Create zip archive " + zipFile.getAbsolutePath(), e);
        }
    }

    public static void compressToZipArchive(OutputStream zipFile, Map<InputStream, String> input) {
        try (
                ZipOutputStream zos = new ZipOutputStream(zipFile)
        ) {
            // create byte buffer
            byte[] buffer = new byte[1024];


            int index = 0;
            for (Map.Entry<InputStream, String> entry : input.entrySet()) {
                String filename = entry.getValue();
                if (filename == null || filename.isEmpty())
                    filename = "file" + index + ".txt";
                index++;
                try (InputStream stream = entry.getKey()) {
                    // begin writing a new ZIP entry, positions the stream to the start of the entry data
                    zos.putNextEntry(new ZipEntry(filename));
                    int length;
                    while ((length = stream.read(buffer)) > 0) {
                        zos.write(buffer, 0, length);

                    }
                    zos.closeEntry();
                    // close the InputStream
                    stream.close();
                } catch (IOException ioe) {
                    LoggerFactory.getLogger(Compress.class).error("Could not Compress " + filename, ioe);
                    throw ioe;
                }
            }

            // close the ZipOutputStream
            zos.close();
        } catch (IOException ioe) {
            LoggerFactory.getLogger(Compress.class).error("Could not Create zip archive", ioe);
        }
    }
}
