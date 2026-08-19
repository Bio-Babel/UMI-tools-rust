//! The `--help` text of each subcommand
//! (Python source: what `optparse` prints for the parser each `main()` builds,
//! as customised by `Utilities.BetterFormatter`).
//!
//! Captured from the live command by
//! `port_reports/umi_tools/data_exports/export_help.py` rather than
//! re-derived. `optparse`'s `IndentedHelpFormatter` has its own wrapping,
//! indent positions and `%default` expansion, and the package's own test suite
//! compares this text byte for byte — so reimplementing the formatter and
//! re-typing some sixty option help strings per subcommand would risk a silent
//! discrepancy for no gain. What is here is exactly what the original prints.
//!
//! `clap` still owns *parsing*; this text is only what `--help` displays. The
//! two are kept in step by `crates/cli`'s tests, which assert that every
//! option named in the help text is one the parser accepts.

/// The `umi_tools extract --help` text.
pub const EXTRACT: &str = r"
extract - Extract UMI from fastq

Usage:

   Single-end:
      umi_tools extract [OPTIONS] -p PATTERN [-I IN_FASTQ[.gz]] [-S OUT_FASTQ[.gz]]

   Paired end:
      umi_tools extract [OPTIONS] -p PATTERN [-I IN_FASTQ[.gz]] [-S OUT_FASTQ[.gz]] --read2-in=IN2_FASTQ[.gz] --read2-out=OUT2_FASTQ[.gz]

   note: If -I/-S are ommited standard in and standard out are used
         for input and output.  To generate a valid BAM file on
         standard out, please redirect log with --log=LOGFILE or
         --log2stderr. Input/Output will be (de)compressed if a
         filename provided to -S/-I/--read2-in/read2-out ends in .gz
         

For full UMI-tools documentation, see https://umi-tools.readthedocs.io/en/latest/

Options:
  --version             show program's version number and exit

  extract-specific options:
    --read2-out=READ2_OUT
                        file to output processed paired read to
    --read2-stdout      Paired reads, send read2 to stdout, discarding read1
    --quality-filter-threshold=QUALITY_FILTER_THRESHOLD
                        Remove reads where any UMI base quality score falls
                        below this threshold
    --quality-filter-mask=QUALITY_FILTER_MASK
                        If a UMI base has a quality below this threshold,
                        replace the base with 'N'
    --quality-encoding=QUALITY_ENCODING
                        Quality score encoding. Choose from 'phred33'[33-77]
                        'phred64' [64-106] or 'solexa' [59-106]
    --error-correct-cell
                        Correct errors in the cell barcode
    --whitelist=WHITELIST
                        A whitelist of accepted cell barcodes
    --blacklist=BLACKLIST
                        A blacklist of rejected cell barcodes
    --subset-reads=READS_SUBSET, --reads-subset=READS_SUBSET
                        Only extract from the first N reads. If N is greater
                        than the number of reads, all reads will be used
    --reconcile-pairs   Allow the presences of reads in read2 input that are
                        not present in read1 input. This allows cell barcode
                        filtering of read1s without considering read2s
    --umi-separator=UMI_SEPARATOR
                        Separator to use to add UMI to the read name. Default:
                        _

  [EXPERIMENTAl] barcode extraction options:
    --either-read       UMI may be on either read (see --either-read-resolve)
                        for options to resolve cases whereUMI is on both reads
    --either-read-resolve=EITHER_READ_RESOLVE
                        How to resolve instances where both reads contain a
                        UMI but using --either-read.Choose from 'discard' or
                        'quality'(use highest quality). default=dicard

  fastq barcode extraction options:
    --extract-method=EXTRACT_METHOD
                        How to extract the umi +/- cell barcodes, Choose from
                        'string' or 'regex'
    -p PATTERN, --bc-pattern=PATTERN
                        Barcode pattern
    --bc-pattern2=PATTERN2
                        Barcode pattern for paired reads
    --3prime            barcode is on 3' end of read.
    --read2-in=READ2_IN
                        file name for read pairs
    --read2-only        Only extract from read2
    --filtered-out=FILTERED_OUT
                        Write out reads not matching regex pattern to this
                        file
    --filtered-out2=FILTERED_OUT2
                        Write out paired reads not matching regex pattern to
                        this file
    --ignore-read-pair-suffixes
                        Ignore '\1' and '\2' read name suffixes

  Input/Output pipe options:
    -I FILE, --stdin=FILE
                        file to read stdin from [default = stdin].
    -L FILE, --log=FILE
                        file with logging information [default = stdout].
    -E FILE, --error=FILE
                        file with error information [default = stderr].
    -S FILE, --stdout=FILE
                        file where output is to go [default = stdout].
    --temp-dir=FILE     Directory for temporary files. If not set, the bash
                        environmental variable TMPDIR is used[default = None].
    --log2stderr        send logging information to stderr [default = False].
    --compresslevel=COMPRESSLEVEL
                        Level of Gzip compression to use. Default (6)
                        matchesGNU gzip rather than python gzip default (which
                        is 9)

  profiling options:
    --timeit=TIMEIT_FILE
                        store timeing information in file [none].
    --timeit-name=TIMEIT_NAME
                        name in timing file for this class of jobs [all].
    --timeit-header     add header for timing information [none].

  common options:
    -v LOGLEVEL, --verbose=LOGLEVEL
                        loglevel [1]. The higher, the more output.
    -h, --help          output short help (command line options only).
    --help-extended     Output full documentation
    --random-seed=RANDOM_SEED
                        random seed to initialize number generator with
                        [none].
";

/// The `umi_tools whitelist --help` text.
pub const WHITELIST: &str = r"
whitelist - Generates a whitelist of accepted cell barcodes

Usage:

   Single-end:
      umi_tools whitelist [OPTIONS] [-I IN_FASTQ[.gz]] [-S OUT_TSV[.gz]]

   Paired end:
      umi_tools whitelist [OPTIONS] [-I IN_FASTQ[.gz]] [-S OUT_TSV[.gz]] --read2-in=IN2_FASTQ[.gz]

   note: If -I/-S are ommited standard in and standard out are used
         for input and output.  Input/Output will be (de)compressed if a
         filename provided to -S/-I/--read2-in ends in .gz


For full UMI-tools documentation, see https://umi-tools.readthedocs.io/en/latest/

Options:
  --version             show program's version number and exit
  --ed-above-threshold=ED_ABOVE_THRESHOLD
                        Detect CBs above the threshold which may be sequence
                        errors from another CB and either 'discard' or
                        'correct'. Default=None (No correction)

  whitelist-specific options:
    --plot-prefix=PLOT_PREFIX
                        Prefix for plots to visualise the automated detection
                        of the number of 'true' cell barcodes
    --subset-reads=SUBSET_READS
                        Use the first N reads to automatically identify the
                        true cell barcodes. If N is greater than the number of
                        reads, all reads will be used. Default is 100,000,000
    --error-correct-threshold=ERROR_CORRECT_THRESHOLD
                        Hamming distance for correction of barcodes to
                        whitelist barcodes. This value will also be used for
                        error detection above the knee if required (--ed-
                        above-threshold)
    --method=METHOD     Use reads or unique umi counts per cell
    --knee-method=KNEE_METHOD
                        Use distance or density methods for detection of knee
    --expect-cells=EXPECT_CELLS
                        Prior expectation on the upper limit on the number of
                        cells sequenced
    --allow-threshold-error
                        Don't select a threshold. Will still output the plots
                        if requested (--plot-prefix)
    --set-cell-number=CELL_NUMBER
                        Specify the number of cell barcodes to accept

  fastq barcode extraction options:
    --extract-method=EXTRACT_METHOD
                        How to extract the umi +/- cell barcodes, Choose from
                        'string' or 'regex'
    -p PATTERN, --bc-pattern=PATTERN
                        Barcode pattern
    --bc-pattern2=PATTERN2
                        Barcode pattern for paired reads
    --3prime            barcode is on 3' end of read.
    --read2-in=READ2_IN
                        file name for read pairs
    --read2-only        Only extract from read2
    --filtered-out=FILTERED_OUT
                        Write out reads not matching regex pattern to this
                        file
    --filtered-out2=FILTERED_OUT2
                        Write out paired reads not matching regex pattern to
                        this file
    --ignore-read-pair-suffixes
                        Ignore '\1' and '\2' read name suffixes

  Input/Output pipe options:
    -I FILE, --stdin=FILE
                        file to read stdin from [default = stdin].
    -L FILE, --log=FILE
                        file with logging information [default = stdout].
    -E FILE, --error=FILE
                        file with error information [default = stderr].
    -S FILE, --stdout=FILE
                        file where output is to go [default = stdout].
    --temp-dir=FILE     Directory for temporary files. If not set, the bash
                        environmental variable TMPDIR is used[default = None].
    --log2stderr        send logging information to stderr [default = False].
    --compresslevel=COMPRESSLEVEL
                        Level of Gzip compression to use. Default (6)
                        matchesGNU gzip rather than python gzip default (which
                        is 9)

  profiling options:
    --timeit=TIMEIT_FILE
                        store timeing information in file [none].
    --timeit-name=TIMEIT_NAME
                        name in timing file for this class of jobs [all].
    --timeit-header     add header for timing information [none].

  common options:
    -v LOGLEVEL, --verbose=LOGLEVEL
                        loglevel [1]. The higher, the more output.
    -h, --help          output short help (command line options only).
    --help-extended     Output full documentation
    --random-seed=RANDOM_SEED
                        random seed to initialize number generator with
                        [none].
";

/// The `umi_tools group --help` text.
pub const GROUP: &str = r"
group - Group reads based on their UMI

Usage: umi_tools group --output-bam [OPTIONS] [--stdin=INFILE.bam] [--stdout=OUTFILE.bam]

       note: If --stdout is ommited, standard out is output. To
             generate a valid BAM file on standard out, please
             redirect log with --log=LOGFILE or --log2stderr 

For full UMI-tools documentation, see https://umi-tools.readthedocs.io/en/latest/

Options:
  --version             show program's version number and exit
  --umi-group-tag=UMI_GROUP_TAG
                        tag for the outputted umi group

  group-specific options:
    --group-out=TSV     Outfile name for file mapping read id to read group
    --output-bam        output a bam file with read groups tagged using the UG
                        tag[default=False]

  Barcode extraction options:
    --extract-umi-method=GET_UMI_METHOD
                        how is the read UMI +/ cell barcode encoded?
                        [default=read_id]
    --umi-separator=UMI_SEP
                        separator between read id and UMI
    --umi-tag=UMI_TAG   tag containing umi
    --umi-tag-split=UMI_TAG_SPLIT
                        split UMI in tag and take the first element
    --umi-tag-delimiter=UMI_TAG_DELIM
                        concatenate UMI in tag separated by delimiter
    --cell-tag=CELL_TAG
                        tag containing cell barcode
    --cell-tag-split=CELL_TAG_SPLIT
                        split cell barcode in tag and take the firstelement
                        for e.g 10X GEM tags
    --cell-tag-delimiter=CELL_TAG_DELIM
                        concatenate cell barcode in tag separated by delimiter

  UMI grouping options:
    --method=METHOD     method to use for umi grouping [default=directional]
    --edit-distance-threshold=THRESHOLD
                        Edit distance theshold at which to join two UMIs when
                        grouping UMIs. [default=1]
    --spliced-is-unique
                        Treat a spliced read as different to an unspliced one
                        [default=False]
    --soft-clip-threshold=SOFT_CLIP_THRESHOLD
                        number of bases clipped from 5' end before read is
                        counted as spliced [default=4]
    --read-length       use read length in addition to position and UMI to
                        identify possible duplicates [default=False]

  single-cell RNA-Seq options:
    --per-gene          Group/Dedup/Count per gene. Must combine with either
                        --gene-tag or --per-contig
    --gene-tag=GENE_TAG
                        Gene is defined by this bam tag [default=none]
    --assigned-status-tag=ASSIGNED_TAG
                        Bam tag describing whether read is assigned to a gene
                        By defualt, this is set as the same tag as --gene-tag
    --skip-tags-regex=SKIP_REGEX
                        Used with --gene-tag. Ignore reads where the gene-tag
                        matches this regex
    --per-contig        group/dedup/count UMIs per contig (field 3 in BAM;
                        RNAME), e.g for transcriptome where contig = gene
    --gene-transcript-map=GENE_TRANSCRIPT_MAP
                        File mapping transcripts to genes (tab separated)
    --per-cell          group/dedup/count per cell

  group/dedup options:
    --buffer-whole-contig
                        Read whole contig before outputting bundles:
                        guarantees that no reads are missed, but increases
                        memory usage
    --multimapping-detection-method=DETECTION_METHOD
                        Some aligners identify multimapping using bam tags.
                        Setting this option to NH, X0 or XT will use these
                        tags when selecting the best read amongst reads with
                        the same position and umi [default=none]

  Input/Output format options:
    --in-format=IN_FORMAT
                        File format of the input file. Format is usually
                        implied from the extension of the filename, but maybe
                        overridden with this option. Default=bam
    --input-options=INPUT_OPTIONS
                        Format string provided to htslib for reading. Mostly
                        useful for CRAM formatted files. See samtools
                        documentation
    -i, --in-sam        [Deprecated] Input file is in sam format
                        [default=False]
    --out-format=OUT_FORMAT
                        File format of the input file. Format is usually
                        implied from the extension of the filename, but maybe
                        overridden with this option. Default=bam
    --output-options=OUTPUT_OPTIONS
                        Format string provided to htslib for writing. Mostly
                        useful for CRAM formatted files. See samtools
                        documentation
    -o, --out-sam       [Deprecated] Output alignments in sam format
                        [default=False]
    --reference-filename=REFERENCE_FILENAME
                        File path or URL to the genome reference to be used
                        when reading or writing CRAM files. By default, when
                        reading a CRAM file, the reference recorded in the
                        input file will be used unless this is specified. When
                        writing, specifying a reference location is required.

  SAM/BAM filtering options:
    --mapping-quality=MAPPING_QUALITY
                        Minimum mapping quality for a read to be retained
                        [default=0]
    --ignore-umi        Ignore UMI and dedup only on position
    --ignore-tlen       Option to dedup paired end reads based solely on
                        read1, whether or not the template length is the same
    --chrom=CHROM       Restrict to one chromosome
    --subset=SUBSET     Use only a fraction of reads, specified by subset
    --paired            paired input BAM. [default=False]
    --no-sort-output    Don't Sort the output

  Group SAM/BAM options:
    --unmapped-reads=UNMAPPED_READS
                        How to handle unmapped reads. Options are 'discard',
                        'use' or 'output' [default=discard]
    --chimeric-pairs=CHIMERIC_PAIRS
                        How to handle chimeric read pairs. Options are
                        'discard', 'use' or 'output' [default=use]
    --unpaired-reads=UNPAIRED_READS
                        How to handle unpaired reads. Options are 'discard',
                        'use' or 'output' [default=use]

  Input/Output pipe options:
    -I FILE, --stdin=FILE
                        file to read stdin from [default = stdin].
    -L FILE, --log=FILE
                        file with logging information [default = stdout].
    -E FILE, --error=FILE
                        file with error information [default = stderr].
    -S FILE, --stdout=FILE
                        file where output is to go [default = stdout].
    --temp-dir=FILE     Directory for temporary files. If not set, the bash
                        environmental variable TMPDIR is used[default = None].
    --log2stderr        send logging information to stderr [default = False].
    --compresslevel=COMPRESSLEVEL
                        Level of Gzip compression to use. Default (6)
                        matchesGNU gzip rather than python gzip default (which
                        is 9)

  profiling options:
    --timeit=TIMEIT_FILE
                        store timeing information in file [none].
    --timeit-name=TIMEIT_NAME
                        name in timing file for this class of jobs [all].
    --timeit-header     add header for timing information [none].

  common options:
    -v LOGLEVEL, --verbose=LOGLEVEL
                        loglevel [1]. The higher, the more output.
    -h, --help          output short help (command line options only).
    --help-extended     Output full documentation
    --random-seed=RANDOM_SEED
                        random seed to initialize number generator with
                        [none].
";

/// The `umi_tools dedup --help` text.
pub const DEDUP: &str = r"
dedup - Deduplicate reads using UMI and mapping coordinates

Usage: umi_tools dedup [OPTIONS] [--stdin=IN_BAM] [--stdout=OUT_BAM]

       note: If --stdout is ommited, standard out is output. To
             generate a valid BAM file on standard out, please
             redirect log with --log=LOGFILE or --log2stderr 

For full UMI-tools documentation, see https://umi-tools.readthedocs.io/en/latest/

Options:
  --version             show program's version number and exit

  dedup-specific options:
    --output-stats=STATS
                        Specify location to output stats

  Barcode extraction options:
    --extract-umi-method=GET_UMI_METHOD
                        how is the read UMI +/ cell barcode encoded?
                        [default=read_id]
    --umi-separator=UMI_SEP
                        separator between read id and UMI
    --umi-tag=UMI_TAG   tag containing umi
    --umi-tag-split=UMI_TAG_SPLIT
                        split UMI in tag and take the first element
    --umi-tag-delimiter=UMI_TAG_DELIM
                        concatenate UMI in tag separated by delimiter
    --cell-tag=CELL_TAG
                        tag containing cell barcode
    --cell-tag-split=CELL_TAG_SPLIT
                        split cell barcode in tag and take the firstelement
                        for e.g 10X GEM tags
    --cell-tag-delimiter=CELL_TAG_DELIM
                        concatenate cell barcode in tag separated by delimiter

  UMI grouping options:
    --method=METHOD     method to use for umi grouping [default=directional]
    --edit-distance-threshold=THRESHOLD
                        Edit distance theshold at which to join two UMIs when
                        grouping UMIs. [default=1]
    --spliced-is-unique
                        Treat a spliced read as different to an unspliced one
                        [default=False]
    --soft-clip-threshold=SOFT_CLIP_THRESHOLD
                        number of bases clipped from 5' end before read is
                        counted as spliced [default=4]
    --read-length       use read length in addition to position and UMI to
                        identify possible duplicates [default=False]

  single-cell RNA-Seq options:
    --per-gene          Group/Dedup/Count per gene. Must combine with either
                        --gene-tag or --per-contig
    --gene-tag=GENE_TAG
                        Gene is defined by this bam tag [default=none]
    --assigned-status-tag=ASSIGNED_TAG
                        Bam tag describing whether read is assigned to a gene
                        By defualt, this is set as the same tag as --gene-tag
    --skip-tags-regex=SKIP_REGEX
                        Used with --gene-tag. Ignore reads where the gene-tag
                        matches this regex
    --per-contig        group/dedup/count UMIs per contig (field 3 in BAM;
                        RNAME), e.g for transcriptome where contig = gene
    --gene-transcript-map=GENE_TRANSCRIPT_MAP
                        File mapping transcripts to genes (tab separated)
    --per-cell          group/dedup/count per cell

  group/dedup options:
    --buffer-whole-contig
                        Read whole contig before outputting bundles:
                        guarantees that no reads are missed, but increases
                        memory usage
    --multimapping-detection-method=DETECTION_METHOD
                        Some aligners identify multimapping using bam tags.
                        Setting this option to NH, X0 or XT will use these
                        tags when selecting the best read amongst reads with
                        the same position and umi [default=none]

  Input/Output format options:
    --in-format=IN_FORMAT
                        File format of the input file. Format is usually
                        implied from the extension of the filename, but maybe
                        overridden with this option. Default=bam
    --input-options=INPUT_OPTIONS
                        Format string provided to htslib for reading. Mostly
                        useful for CRAM formatted files. See samtools
                        documentation
    -i, --in-sam        [Deprecated] Input file is in sam format
                        [default=False]
    --out-format=OUT_FORMAT
                        File format of the input file. Format is usually
                        implied from the extension of the filename, but maybe
                        overridden with this option. Default=bam
    --output-options=OUTPUT_OPTIONS
                        Format string provided to htslib for writing. Mostly
                        useful for CRAM formatted files. See samtools
                        documentation
    -o, --out-sam       [Deprecated] Output alignments in sam format
                        [default=False]
    --reference-filename=REFERENCE_FILENAME
                        File path or URL to the genome reference to be used
                        when reading or writing CRAM files. By default, when
                        reading a CRAM file, the reference recorded in the
                        input file will be used unless this is specified. When
                        writing, specifying a reference location is required.

  SAM/BAM filtering options:
    --mapping-quality=MAPPING_QUALITY
                        Minimum mapping quality for a read to be retained
                        [default=0]
    --ignore-umi        Ignore UMI and dedup only on position
    --ignore-tlen       Option to dedup paired end reads based solely on
                        read1, whether or not the template length is the same
    --chrom=CHROM       Restrict to one chromosome
    --subset=SUBSET     Use only a fraction of reads, specified by subset
    --paired            paired input BAM. [default=False]
    --no-sort-output    Don't Sort the output

  Dedup and Count SAM/BAM options:
    --unmapped-reads=UNMAPPED_READS
                        How to handle unmapped reads. Options are 'discard' or
                        'use' [default=discard]
    --chimeric-pairs=CHIMERIC_PAIRS
                        How to handle chimeric read pairs. Options are
                        'discard' or 'use'  [default=use]
    --unpaired-reads=UNPAIRED_READS
                        How to handle unpaired reads. Options are 'discard'or
                        'use' [default=use]

  Input/Output pipe options:
    -I FILE, --stdin=FILE
                        file to read stdin from [default = stdin].
    -L FILE, --log=FILE
                        file with logging information [default = stdout].
    -E FILE, --error=FILE
                        file with error information [default = stderr].
    -S FILE, --stdout=FILE
                        file where output is to go [default = stdout].
    --temp-dir=FILE     Directory for temporary files. If not set, the bash
                        environmental variable TMPDIR is used[default = None].
    --log2stderr        send logging information to stderr [default = False].
    --compresslevel=COMPRESSLEVEL
                        Level of Gzip compression to use. Default (6)
                        matchesGNU gzip rather than python gzip default (which
                        is 9)

  profiling options:
    --timeit=TIMEIT_FILE
                        store timeing information in file [none].
    --timeit-name=TIMEIT_NAME
                        name in timing file for this class of jobs [all].
    --timeit-header     add header for timing information [none].

  common options:
    -v LOGLEVEL, --verbose=LOGLEVEL
                        loglevel [1]. The higher, the more output.
    -h, --help          output short help (command line options only).
    --help-extended     Output full documentation
    --random-seed=RANDOM_SEED
                        random seed to initialize number generator with
                        [none].
";

/// The `umi_tools count --help` text.
pub const COUNT: &str = r"
count - Count reads per-gene using UMI and mapping coordinates

Usage: umi_tools count [OPTIONS] --stdin=IN_BAM [--stdout=OUT_BAM]

       note: If --stdout is ommited, standard out is output. To
             generate a valid BAM file on standard out, please
             redirect log with --log=LOGFILE or --log2stderr 

For full UMI-tools documentation, see https://umi-tools.readthedocs.io/en/latest/

Options:
  --version             show program's version number and exit
  --wide-format-cell-counts
                        output the cell counts in a wide format (rows=genes,
                        columns=cells)

  count-specific options:

  Barcode extraction options:
    --extract-umi-method=GET_UMI_METHOD
                        how is the read UMI +/ cell barcode encoded?
                        [default=read_id]
    --umi-separator=UMI_SEP
                        separator between read id and UMI
    --umi-tag=UMI_TAG   tag containing umi
    --umi-tag-split=UMI_TAG_SPLIT
                        split UMI in tag and take the first element
    --umi-tag-delimiter=UMI_TAG_DELIM
                        concatenate UMI in tag separated by delimiter
    --cell-tag=CELL_TAG
                        tag containing cell barcode
    --cell-tag-split=CELL_TAG_SPLIT
                        split cell barcode in tag and take the firstelement
                        for e.g 10X GEM tags
    --cell-tag-delimiter=CELL_TAG_DELIM
                        concatenate cell barcode in tag separated by delimiter

  UMI grouping options:
    --method=METHOD     method to use for umi grouping [default=directional]
    --edit-distance-threshold=THRESHOLD
                        Edit distance theshold at which to join two UMIs when
                        grouping UMIs. [default=1]
    --spliced-is-unique
                        Treat a spliced read as different to an unspliced one
                        [default=False]
    --soft-clip-threshold=SOFT_CLIP_THRESHOLD
                        number of bases clipped from 5' end before read is
                        counted as spliced [default=4]
    --read-length       use read length in addition to position and UMI to
                        identify possible duplicates [default=False]

  single-cell RNA-Seq options:
    --per-gene          Group/Dedup/Count per gene. Must combine with either
                        --gene-tag or --per-contig
    --gene-tag=GENE_TAG
                        Gene is defined by this bam tag [default=none]
    --assigned-status-tag=ASSIGNED_TAG
                        Bam tag describing whether read is assigned to a gene
                        By defualt, this is set as the same tag as --gene-tag
    --skip-tags-regex=SKIP_REGEX
                        Used with --gene-tag. Ignore reads where the gene-tag
                        matches this regex
    --per-contig        group/dedup/count UMIs per contig (field 3 in BAM;
                        RNAME), e.g for transcriptome where contig = gene
    --gene-transcript-map=GENE_TRANSCRIPT_MAP
                        File mapping transcripts to genes (tab separated)
    --per-cell          group/dedup/count per cell

  Input/Output format options:
    --in-format=IN_FORMAT
                        File format of the input file. Format is usually
                        implied from the extension of the filename, but maybe
                        overridden with this option. Default=bam
    --input-options=INPUT_OPTIONS
                        Format string provided to htslib for reading. Mostly
                        useful for CRAM formatted files. See samtools
                        documentation
    -i, --in-sam        [Deprecated] Input file is in sam format
                        [default=False]
    --reference-filename=REFERENCE_FILENAME
                        File path or URL to the genome reference to be used
                        when reading or writing CRAM files. By default, when
                        reading a CRAM file, the reference recorded in the
                        input file will be used unless this is specified. When
                        writing, specifying a reference location is required.

  SAM/BAM filtering options:
    --mapping-quality=MAPPING_QUALITY
                        Minimum mapping quality for a read to be retained
                        [default=0]
    --ignore-umi        Ignore UMI and dedup only on position
    --ignore-tlen       Option to dedup paired end reads based solely on
                        read1, whether or not the template length is the same
    --chrom=CHROM       Restrict to one chromosome
    --subset=SUBSET     Use only a fraction of reads, specified by subset
    --paired            paired input BAM. [default=False]
    --no-sort-output    Don't Sort the output

  Dedup and Count SAM/BAM options:
    --unmapped-reads=UNMAPPED_READS
                        How to handle unmapped reads. Options are 'discard' or
                        'use' [default=discard]
    --chimeric-pairs=CHIMERIC_PAIRS
                        How to handle chimeric read pairs. Options are
                        'discard' or 'use'  [default=use]
    --unpaired-reads=UNPAIRED_READS
                        How to handle unpaired reads. Options are 'discard'or
                        'use' [default=use]

  Input/Output pipe options:
    -I FILE, --stdin=FILE
                        file to read stdin from [default = stdin].
    -L FILE, --log=FILE
                        file with logging information [default = stdout].
    -E FILE, --error=FILE
                        file with error information [default = stderr].
    -S FILE, --stdout=FILE
                        file where output is to go [default = stdout].
    --temp-dir=FILE     Directory for temporary files. If not set, the bash
                        environmental variable TMPDIR is used[default = None].
    --log2stderr        send logging information to stderr [default = False].
    --compresslevel=COMPRESSLEVEL
                        Level of Gzip compression to use. Default (6)
                        matchesGNU gzip rather than python gzip default (which
                        is 9)

  profiling options:
    --timeit=TIMEIT_FILE
                        store timeing information in file [none].
    --timeit-name=TIMEIT_NAME
                        name in timing file for this class of jobs [all].
    --timeit-header     add header for timing information [none].

  common options:
    -v LOGLEVEL, --verbose=LOGLEVEL
                        loglevel [1]. The higher, the more output.
    -h, --help          output short help (command line options only).
    --help-extended     Output full documentation
    --random-seed=RANDOM_SEED
                        random seed to initialize number generator with
                        [none].
";

/// The `umi_tools count\_tab --help` text.
pub const COUNT_TAB: &str = r"
count_tab - Count reads per gene from flatfile using UMIs

Usage: umi_tools count_tab [OPTIONS] [--stdin=IN_TSV[.gz]] [--stdout=OUT_TSV[.gz]]

       note: If --stdin/--stdout are ommited standard in and standard
             out are used for input and output. Input/Output will be
             (de)compressed if a filename provided to --stdin/--stdout
             ends in .gz 

For full UMI-tools documentation, see https://umi-tools.readthedocs.io/en/latest/

Options:
  --version             show program's version number and exit

  count_tab-specific options:
    --barcode-separator=BC_SEP
                        separator between read id and UMI  and (optionally)
                        the cell barcode
    --per-cell          Readname includes cell barcode as well as UMI in
                        format: read[sep]UMI[sep]CB

  UMI grouping options:
    --method=METHOD     method to use for umi grouping [default=directional]
    --edit-distance-threshold=THRESHOLD
                        Edit distance theshold at which to join two UMIs when
                        grouping UMIs. [default=1]
    --spliced-is-unique
                        Treat a spliced read as different to an unspliced one
                        [default=False]
    --soft-clip-threshold=SOFT_CLIP_THRESHOLD
                        number of bases clipped from 5' end before read is
                        counted as spliced [default=4]
    --read-length       use read length in addition to position and UMI to
                        identify possible duplicates [default=False]

  Input/Output format options:
    --in-format=IN_FORMAT
                        File format of the input file. Format is usually
                        implied from the extension of the filename, but maybe
                        overridden with this option. Default=bam
    --input-options=INPUT_OPTIONS
                        Format string provided to htslib for reading. Mostly
                        useful for CRAM formatted files. See samtools
                        documentation
    -i, --in-sam        [Deprecated] Input file is in sam format
                        [default=False]
    --reference-filename=REFERENCE_FILENAME
                        File path or URL to the genome reference to be used
                        when reading or writing CRAM files. By default, when
                        reading a CRAM file, the reference recorded in the
                        input file will be used unless this is specified. When
                        writing, specifying a reference location is required.

  Input/Output pipe options:
    -I FILE, --stdin=FILE
                        file to read stdin from [default = stdin].
    -L FILE, --log=FILE
                        file with logging information [default = stdout].
    -E FILE, --error=FILE
                        file with error information [default = stderr].
    -S FILE, --stdout=FILE
                        file where output is to go [default = stdout].
    --temp-dir=FILE     Directory for temporary files. If not set, the bash
                        environmental variable TMPDIR is used[default = None].
    --log2stderr        send logging information to stderr [default = False].
    --compresslevel=COMPRESSLEVEL
                        Level of Gzip compression to use. Default (6)
                        matchesGNU gzip rather than python gzip default (which
                        is 9)

  profiling options:
    --timeit=TIMEIT_FILE
                        store timeing information in file [none].
    --timeit-name=TIMEIT_NAME
                        name in timing file for this class of jobs [all].
    --timeit-header     add header for timing information [none].

  common options:
    -v LOGLEVEL, --verbose=LOGLEVEL
                        loglevel [1]. The higher, the more output.
    -h, --help          output short help (command line options only).
    --help-extended     Output full documentation
    --random-seed=RANDOM_SEED
                        random seed to initialize number generator with
                        [none].
";

/// The `umi_tools prepare\_for\_em --help` text.
pub const PREPARE_FOR_EM: &str = r"
prepare_for_em - make output from dedup or group compatible with EM tools

Usage: umi_tools prepare_for_em [OPTIONS] [--stdin=IN_BAM] [--stdout=OUT_BAM]

       note: If --stdout is ommited, standard out is output. To
             generate a valid BAM file on standard out, please
             redirect log with --log=LOGFILE or --log2stderr 

For full UMI-tools documentation, see https://umi-tools.readthedocs.io/en/latest/

Options:
  --version             show program's version number and exit

  RSEM preparation specific options:
    --tags=TAGS         Comma-seperated list of tags to transfer from read1 to
                        read2
    --sam               input and output SAM rather than BAM

  Input/Output format options:
    --in-format=IN_FORMAT
                        File format of the input file. Format is usually
                        implied from the extension of the filename, but maybe
                        overridden with this option. Default=bam
    --input-options=INPUT_OPTIONS
                        Format string provided to htslib for reading. Mostly
                        useful for CRAM formatted files. See samtools
                        documentation
    -i, --in-sam        [Deprecated] Input file is in sam format
                        [default=False]
    --out-format=OUT_FORMAT
                        File format of the input file. Format is usually
                        implied from the extension of the filename, but maybe
                        overridden with this option. Default=bam
    --output-options=OUTPUT_OPTIONS
                        Format string provided to htslib for writing. Mostly
                        useful for CRAM formatted files. See samtools
                        documentation
    -o, --out-sam       [Deprecated] Output alignments in sam format
                        [default=False]
    --reference-filename=REFERENCE_FILENAME
                        File path or URL to the genome reference to be used
                        when reading or writing CRAM files. By default, when
                        reading a CRAM file, the reference recorded in the
                        input file will be used unless this is specified. When
                        writing, specifying a reference location is required.

  Input/Output pipe options:
    -I FILE, --stdin=FILE
                        file to read stdin from [default = stdin].
    -L FILE, --log=FILE
                        file with logging information [default = stdout].
    -E FILE, --error=FILE
                        file with error information [default = stderr].
    -S FILE, --stdout=FILE
                        file where output is to go [default = stdout].
    --temp-dir=FILE     Directory for temporary files. If not set, the bash
                        environmental variable TMPDIR is used[default = None].
    --log2stderr        send logging information to stderr [default = False].
    --compresslevel=COMPRESSLEVEL
                        Level of Gzip compression to use. Default (6)
                        matchesGNU gzip rather than python gzip default (which
                        is 9)

  profiling options:
    --timeit=TIMEIT_FILE
                        store timeing information in file [none].
    --timeit-name=TIMEIT_NAME
                        name in timing file for this class of jobs [all].
    --timeit-header     add header for timing information [none].

  common options:
    -v LOGLEVEL, --verbose=LOGLEVEL
                        loglevel [1]. The higher, the more output.
    -h, --help          output short help (command line options only).
    --help-extended     Output full documentation
    --random-seed=RANDOM_SEED
                        random seed to initialize number generator with
                        [none].
";

/// The help text for a subcommand, by the name the dispatcher uses.
#[must_use]
pub fn for_command(name: &str) -> Option<&'static str> {
    Some(match name {
        "extract" => EXTRACT,
        "whitelist" => WHITELIST,
        "group" => GROUP,
        "dedup" => DEDUP,
        "count" => COUNT,
        "count_tab" => COUNT_TAB,
        "prepare-for-em" | "prepare_for_em" => PREPARE_FOR_EM,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_subcommand_has_help() {
        for name in [
            "extract",
            "whitelist",
            "group",
            "dedup",
            "count",
            "count_tab",
            "prepare_for_em",
            "prepare-for-em",
        ] {
            let text = for_command(name).unwrap_or_else(|| panic!("no help for {name}"));
            assert!(!text.trim().is_empty(), "{name} help is empty");
        }
        assert!(for_command("nonsense").is_none());
    }

    /// Each help text must describe the command it belongs to, which catches a
    /// copy-paste in the generated table.
    #[test]
    fn help_texts_are_not_swapped() {
        assert!(EXTRACT.contains("extract - "));
        assert!(WHITELIST.contains("whitelist - "));
        assert!(GROUP.contains("group - "));
        assert!(DEDUP.contains("dedup - "));
        assert!(COUNT.contains("count - "));
        assert!(COUNT_TAB.contains("count_tab - "));
        assert!(PREPARE_FOR_EM.contains("prepare_for_em"));
    }

    /// The parser and the help text must agree on which options exist: a
    /// flag documented here but unknown to `clap` would be accepted in the
    /// help and rejected on use.
    #[test]
    fn documented_options_are_accepted_by_the_parser() {
        use crate::options::{build_parser, OptionGroups};

        // The full option surface: every group, so any documented flag is
        // covered by some subcommand's parser.
        let command = build_parser(
            "all",
            "every option group",
            OptionGroups {
                extract: true,
                group_dedup: true,
                sam: true,
                dedup_count_sam: true,
                group_sam: true,
                umi_grouping: true,
                sbcram: true,
                exclude_output_format: false,
            },
        );
        // Aliases count as accepted: `--buffer-whole-contig` is an alias of
        // `--whole-contig`, exactly as the original declares two options
        // writing one destination.
        let known: std::collections::BTreeSet<String> = command
            .get_arguments()
            .flat_map(|a| {
                a.get_long()
                    .into_iter()
                    .chain(a.get_all_aliases().into_iter().flatten())
                    .map(str::to_string)
            })
            .collect();

        // Options declared by an individual subcommand rather than a shared
        // group; each is added in that subcommand's own module.
        let subcommand_local: std::collections::BTreeSet<&str> = [
            "whitelist",
            "blacklist",
            "error-correct-cell",
            "read2-out",
            "read2-stdout",
            "read2-only",
            "retain-umi",
            "subset-reads",
            "reads-subset",
            "reconcile-pairs",
            "umi-separator",
            "either-read",
            "either-read-resolve",
            "umi-correct-log",
            "plot-prefix",
            "error-correct-threshold",
            "knee-method",
            "expect-cells",
            "allow-threshold-error",
            "set-cell-number",
            "ed-above-threshold",
            "group-out",
            "output-bam",
            "umi-group-tag",
            "output-stats",
            "wide-format-cell-counts",
            "barcode-separator",
            "tags",
            "sam",
            "help",
            "version",
            "short-help",
            "no-usage",
            "test",
            "timeit",
            "timeit-name",
            "timeit-header",
            "filter-umi",
            "umi-whitelist",
            "umi-whitelist-paired",
            "correct-umi-threshold",
        ]
        .into_iter()
        .collect();

        let mut missing: Vec<String> = Vec::new();
        for text in [
            EXTRACT,
            WHITELIST,
            GROUP,
            DEDUP,
            COUNT,
            COUNT_TAB,
            PREPARE_FOR_EM,
        ] {
            for token in text.split_whitespace() {
                let Some(flag) = token.strip_prefix("--") else {
                    continue;
                };
                let flag = flag
                    .trim_end_matches([',', '.', ')', ']', '`', '\'', '"', ':', ';'])
                    .split('=')
                    .next()
                    .unwrap_or("");
                if flag.is_empty()
                    || !flag
                        .chars()
                        .all(|c| c.is_ascii_lowercase() || c == '-' || c.is_ascii_digit())
                {
                    continue;
                }
                if !known.contains(flag) && !subcommand_local.contains(flag) {
                    missing.push(flag.to_string());
                }
            }
        }
        missing.sort();
        missing.dedup();
        assert!(
            missing.is_empty(),
            "help documents options the parser does not accept: {missing:?}"
        );
    }
}
