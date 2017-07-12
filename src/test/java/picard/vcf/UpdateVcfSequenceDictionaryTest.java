/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.vcf;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.io.*;


/**
 * @author George Grant
 */
public class UpdateVcfSequenceDictionaryTest {
    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");
    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("UpdateVcfSequenceDictionaryTest", null);
    private static final File UNZIPPED_FILE = new File(OUTPUT_DATA_PATH + "out.vcf");

    private static final String STD_OUT_NAME = "/dev/stdout";

    @AfterClass
    public void teardown() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @DataProvider(name = "OutputFiles")
    public static Object[][] outputFies() {

        return new Object[][] {
                {OUTPUT_DATA_PATH + "updateVcfSequenceDictionaryTest-delete-me"  + VcfUtils.COMPRESSED_VCF_ENDING},
                {OUTPUT_DATA_PATH + "updateVcfSequenceDictionaryTest-delete-me" + VcfUtils.UNCOMPRESSED_VCF_ENDING},
                {STD_OUT_NAME}
        };
    }

    /**
     * Utility for unzipping a zipped file's contents into a human readable, unzipped file
     *
     * @param zippedFile    input zipped file
     * @param unzippedFile  unzipped file
     * @throws IOException
     */
    private void unzipFile(final File zippedFile, final File unzippedFile) throws IOException {
        final InputStream gzInputStream = IOUtil.openFileForReading(zippedFile);
        final FileOutputStream fileOutputStream = new FileOutputStream(unzippedFile.getAbsolutePath());
        final byte[] buffer = new byte[1024];
        int len;
        while ((len = gzInputStream.read(buffer)) > 0) {
            fileOutputStream.write(buffer, 0, len);
        }
        gzInputStream.close();
        fileOutputStream.close();
    }

    @Test(dataProvider = "OutputFiles")
    public void testUpdateVcfSequenceDictionary(final String outputFileName) throws IOException, NoSuchFieldException, IllegalAccessException {
        final File inputFile = new File(TEST_DATA_PATH, "vcfFormatTest.vcf");
        // vcfFormatTest.bad_dict.vcf is a vcf with two (2) ##contig lines deleted
        final File samSequenceDictionaryVcf = new File(TEST_DATA_PATH, "vcfFormatTest.bad_dict.vcf");
        final File outputFile = new File(outputFileName);
        outputFile.deleteOnExit();
        UNZIPPED_FILE.deleteOnExit();

        final UpdateVcfSequenceDictionary updateVcfSequenceDictionary = new UpdateVcfSequenceDictionary();

        updateVcfSequenceDictionary.INPUT = inputFile;
        updateVcfSequenceDictionary.OUTPUT = outputFile;
        updateVcfSequenceDictionary.SEQUENCE_DICTIONARY = samSequenceDictionaryVcf;

        Assert.assertEquals(updateVcfSequenceDictionary.instanceMain(new String[0]), 0);

        // Unzip vcf.gz file
        if (outputFileName.endsWith(VcfUtils.COMPRESSED_VCF_ENDING)) {
            unzipFile(outputFile, UNZIPPED_FILE);
        }

        // Check that the output is equal to the input file if the output is not going to stdout
        if (outputFileName != STD_OUT_NAME) {
            IOUtil.assertFilesEqual(samSequenceDictionaryVcf,
                    outputFileName.endsWith(VcfUtils.COMPRESSED_VCF_ENDING) ? UNZIPPED_FILE : outputFile);
        }
    }
}
