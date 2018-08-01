package ca.on.oicr.pde.deciders;

import com.google.common.collect.Sets;
import java.util.*;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.module.ReturnValue.ExitStatus;
import net.sourceforge.seqware.common.util.Log;

/**
 *
 * @author abenawra@oicr.on.ca
 */
public class MiXCRDecider extends OicrDecider {

    private String produce_transcriptome_bam = "true";
    private String numOfThreads = "6";
    private String mixcrMemory = "24";
    private String queue = "";

    private static final String OICR = "OICR";
    private static final String ILLUMINA = "Illumina"; //If we don't have this passed as parameter, we assume Illumina
    private Set<String> allowedTemplateTypes;

    private String input_read1_fastq;
    private String input_read2_fastq;
    private String external_name;
    private ReadGroupData readGroupDataForWorkflowRun;

    public MiXCRDecider() {
        super();
        parser.accepts("ini-file", "Optional: the location of the INI file.").withRequiredArg();

        //star aln
        parser.accepts("mixcr-mem", "Optional: MiXCR allocated memory Gb, default is 24.").withRequiredArg();

     
        parser.accepts("template-type", "Optional: limit the run to only specified template type(s) (comma separated list).").withRequiredArg();

    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(Arrays.asList("chemical/seq-na-fastq-gzip"));
        this.setHeadersToGroupBy(Arrays.asList(FindAllTheFiles.Header.IUS_SWA));

        //allows anything defined on the command line to override the defaults here.
       //mixcr
        if (this.options.has("mixcr-mem")) {
            this.mixcrMemory = options.valueOf("mixcr-mem").toString();
        }
        if (this.options.has("template-type")) {
            String templateTypeArg = this.options.valueOf("template-type").toString();
            allowedTemplateTypes = Sets.newHashSet(templateTypeArg.split(","));
        }

        ReturnValue val = super.init();

        return val;
    }

    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        this.input_read1_fastq = null;
        this.input_read2_fastq = null;

        String[] filePaths = commaSeparatedFilePaths.split(",");
        if (filePaths.length != 2) {
            Log.error("This Decider supports only cases where we have only 2 files per lane, WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }

        String[] fqFilesArray = commaSeparatedFilePaths.split(",");
        for (String file : fqFilesArray) {
            int mate = idMate(file);
            switch (mate) {
                case 1:
                    if (this.input_read1_fastq != null) {
                        Log.error("More than one file found for read 1: " + input_read1_fastq + ", " + file);
                        return new ReturnValue(ExitStatus.INVALIDFILE);
                    }
                    this.input_read1_fastq = file;
                    break;
                case 2:
                    if (this.input_read2_fastq != null) {
                        Log.error("More than one file found for read 2: " + input_read2_fastq + ", " + file);
                        return new ReturnValue(ExitStatus.INVALIDFILE);
                    }
                    this.input_read2_fastq = file;
                    break;
                default:
                    Log.error("Cannot identify " + file + " end (read 1 or 2)");
                    return new ReturnValue(ExitStatus.INVALIDFILE);
            }
        }

        if (input_read1_fastq == null || input_read2_fastq == null) {
            Log.error("The Decider was not able to find both R1 and R2 fastq files for paired sequencing alignment, WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }

        readGroupDataForWorkflowRun = new ReadGroupData(files.get(input_read1_fastq), files.get(input_read2_fastq));

        return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);

        if (allowedTemplateTypes != null) {
            String currentTemplateType = returnValue.getAttribute(FindAllTheFiles.Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
            if (!allowedTemplateTypes.contains(currentTemplateType)) {
                return false;
            }
        }

        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        Log.debug("INI FILE:" + commaSeparatedFilePaths);

        Map<String, String> iniFileMap = super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        iniFileMap.put("input_file_1", this.input_read1_fastq);
        iniFileMap.put("input_file_2", this.input_read2_fastq);

        iniFileMap.put("mixcr_mem", this.mixcrMemory);
        iniFileMap.put("external_name", this.external_name);

        return iniFileMap;
    }

    public static void main(String args[]) {

        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(MiXCRDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }

}