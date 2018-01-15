/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package picard.sam;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.Md5CalculatingOutputStream;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.ReferenceProgramGroup;
import picard.util.SequenceDictionaryUtils;

import java.io.*;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

/**
 * Create a SAM/BAM file from a fasta containing reference sequence. The output SAM file contains a header but no
 * SAMRecords, and the header contains only sequence records.
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = CreateSequenceDictionary.USAGE_SUMMARY + CreateSequenceDictionary.USAGE_DETAILS,
        oneLineSummary = CreateSequenceDictionary.USAGE_SUMMARY,
        programGroup = ReferenceProgramGroup.class
)
public class CreateSequenceDictionary extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Creates a sequence dictionary for a reference sequence.  ";
    static final String USAGE_DETAILS = "This tool creates a sequence dictionary file (with \".dict\" extension) from a reference " +
            "sequence provided in FASTA format, which is required by many processing and analysis tools. The output file contains a " +
            "header but no SAMRecords, and the header contains only sequence records." +
            "<br /><br />" +
            "The reference sequence can be gzipped (both .fasta and .fasta.gz are supported)." +
            "" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CreateSequenceDictionary \\ <br />" +
            "      R=reference.fasta \\ <br />" +
            "      O=reference.dict" +
            "" +
            "</pre>" +
            "<hr />";
    // The following attributes define the command-line arguments

    private static final Log logger = Log.getInstance(CreateSequenceDictionary.class);

    @Argument(doc = "Output SAM file containing only the sequence dictionary. By default it will use the base name of the input reference with the .dict extension",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true)
    public File OUTPUT;

    @Argument(shortName = "AS", doc = "Put into AS field of sequence dictionary entry if supplied", optional = true)
    public String GENOME_ASSEMBLY;

    @Argument(shortName = "UR", doc = "Put into UR field of sequence dictionary entry.  If not supplied, input reference file is used",
            optional = true)
    public String URI;

    @Argument(shortName = "SP", doc = "Put into SP field of sequence dictionary entry", optional = true)
    public String SPECIES;

    @Argument(doc = "Make sequence name the first word from the > line in the fasta file.  " +
            "By default the entire contents of the > line is used, excluding leading and trailing whitespace.")
    public boolean TRUNCATE_NAMES_AT_WHITESPACE = true;

    @Argument(doc = "Stop after writing this many sequences.  For testing.")
    public int NUM_SEQUENCES = Integer.MAX_VALUE;

    private final MessageDigest md5;

    public CreateSequenceDictionary() {
        try {
            md5 = MessageDigest.getInstance("MD5");
        } catch (NoSuchAlgorithmException e) {
            throw new PicardException("MD5 algorithm not found", e);
        }
    }

    /**
     * Read all the sequences from the given reference file, and convert into SAMSequenceRecords
     *
     * @param referenceFile fasta or fasta.gz
     * @return SAMSequenceRecords containing info from the fasta, plus from cmd-line arguments.
     *
     * @deprecated 12/9/16
     */
    @Deprecated
    public SAMSequenceDictionary makeSequenceDictionary(final File referenceFile) {
        final ReferenceSequenceFile refSeqFile =
                ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile, TRUNCATE_NAMES_AT_WHITESPACE);
        ReferenceSequence refSeq;
        final List<SAMSequenceRecord> ret = new ArrayList<>();
        final Set<String> sequenceNames = new HashSet<>();
        for (int numSequences = 0; numSequences < NUM_SEQUENCES && (refSeq = refSeqFile.nextSequence()) != null; ++numSequences) {
            if (sequenceNames.contains(refSeq.getName())) {
                throw new PicardException("Sequence name appears more than once in reference: " + refSeq.getName());
            }
            sequenceNames.add(refSeq.getName());
            ret.add(makeSequenceRecord(refSeq));
        }
        return new SAMSequenceDictionary(ret);
    }

    /**
     * Use reference filename to create URI to go into header if URI was not passed on cmd line.
     */
    protected String[] customCommandLineValidation() {
        if (URI == null) {
            URI = "file:" + referenceSequence.getReferenceFile().getAbsolutePath();
        }
        if (OUTPUT == null) {
            OUTPUT = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(referenceSequence.getReferenceFile());
            logger.info("Output dictionary will be written in ", OUTPUT);
        }
        return super.customCommandLineValidation();
    }

    // return a custom argument collection because this tool uses the argument name
    // "REFERENCE" instead of "REFERENCE_SEQUENCE"
    @Override
    protected ReferenceArgumentCollection makeReferenceArgumentCollection() {
        return new CreateSeqDictReferenceArgumentCollection();
    }

    public static class CreateSeqDictReferenceArgumentCollection implements ReferenceArgumentCollection {
        @Argument(doc = "Input reference fasta or fasta.gz", shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME)
        public File REFERENCE;

        @Override
        public File getReferenceFile() {
            return REFERENCE;
        }
    }

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    protected int doWork() {
        if (OUTPUT.exists()) {
            throw new PicardException(OUTPUT.getAbsolutePath() +
                    " already exists.  Delete this file and try again, or specify a different output file.");
        }

        try (BufferedWriter writer = makeWriter()) {
            SequenceDictionaryUtils.encodeDictionary(writer, new SAMSequenceRecordIterator(REFERENCE_SEQUENCE, TRUNCATE_NAMES_AT_WHITESPACE));
        } catch (FileNotFoundException e) {
            throw new PicardException("File " + OUTPUT.getAbsolutePath() + " not found");
        } catch (IOException e) {
            throw new PicardException("Can't write to or close output file " + OUTPUT.getAbsolutePath());
        } catch (IllegalArgumentException e) {
            OUTPUT.delete();
            throw new PicardException(e.getMessage());
        }
        return 0;
    }

    private BufferedWriter makeWriter() throws FileNotFoundException {
        return new BufferedWriter(
                new AsciiWriter(this.CREATE_MD5_FILE ?
                        new Md5CalculatingOutputStream(
                                new FileOutputStream(OUTPUT, false),
                                new File(OUTPUT.getAbsolutePath() + ".md5")
                        )
                        : new FileOutputStream(OUTPUT)
                )
        );
    }

    /**
     * Create one SAMSequenceRecord from a single fasta sequence
     */
    private SAMSequenceRecord makeSequenceRecord(final ReferenceSequence refSeq) {
        final SAMSequenceRecord ret = new SAMSequenceRecord(refSeq.getName(), refSeq.length());

        // Compute MD5 of upcased bases
        final byte[] bases = refSeq.getBases();
        for (int i = 0; i < bases.length; ++i) {
            bases[i] = StringUtil.toUpperCase(bases[i]);
        }

        ret.setAttribute(SAMSequenceRecord.MD5_TAG, md5Hash(bases));
        if (GENOME_ASSEMBLY != null) {
            ret.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, GENOME_ASSEMBLY);
        }
        ret.setAttribute(SAMSequenceRecord.URI_TAG, URI);
        if (SPECIES != null) {
            ret.setAttribute(SAMSequenceRecord.SPECIES_TAG, SPECIES);
        }
        return ret;
    }

    private String md5Hash(final byte[] bytes) {
        md5.reset();
        md5.update(bytes);
        String s = new BigInteger(1, md5.digest()).toString(16);
        if (s.length() != 32) {
            final String zeros = "00000000000000000000000000000000";
            s = zeros.substring(0, 32 - s.length()) + s;
        }
        return s;
    }

    private class SAMSequenceRecordIterator implements Iterator<SAMSequenceRecord> {

        private final ReferenceSequenceFile refSeqFile;
        private SAMSequenceRecord nextRecord;

        SAMSequenceRecordIterator(final File reference, boolean truncateAtWhiteSpace) {
            refSeqFile = ReferenceSequenceFileFactory.
                    getReferenceSequenceFile(reference, truncateAtWhiteSpace);
            getNext();
        }

        private void getNext() {
            nextRecord = Optional.ofNullable(refSeqFile.nextSequence())
                    .map(CreateSequenceDictionary.this::makeSequenceRecord)
                    .orElse(null);
        }

        @Override
        public SAMSequenceRecord next() {
            final SAMSequenceRecord tempNext = nextRecord;
            getNext();
            return tempNext;
        }

        @Override
        public boolean hasNext() {
            return nextRecord != null;
        }
    }
}
