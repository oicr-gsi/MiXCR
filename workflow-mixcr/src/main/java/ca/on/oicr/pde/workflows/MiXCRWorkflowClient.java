package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

/**
 * <p>
 * For more information on developing workflows, see the documentation at
 * <a href="http://seqware.github.io/docs/6-pipeline/java-workflows/">SeqWare
 * Java Workflows</a>.</p>
 *
 * Quick reference for the order of methods called: 1. setupDirectory 2.
 * setupFiles 3. setupWorkflow 4. setupvironment 5. buildWorkflow
 *
 * See the SeqWare API for
 * <a href="http://seqware.github.io/javadoc/stable/apidocs/net/sourceforge/seqware/pipeline/workflowV2/AbstractWorkflowDataModel.html#setupDirectory%28%29">AbstractWorkflowDataModel</a>
 * for more information.
 */
public class MiXCRWorkflowClient extends OicrWorkflow {

    //dir
    private String dataDir, tmpDir;
    private String outDir;

    // Input Data
    private String read1Fastq;
    private String read2Fastq;
    private String outputFilenamePrefix;

    //varscan intermediate file names
    private String alignVdjcaFile;
    private String alignRescued1File;
    private String alignRescued2File;
    private String alignRescued2ExtendFile;
    private String cloneClnsFile;
    private String cloneDetTxtFile;

    private String exports;

    //Tools
    private String mixcr;
    private String javaHome;

    //Memory allocation
    private Integer mixcrMem;

    //path to bin
    private String bin;

    private boolean manualOutput;
    private static final Logger logger = Logger.getLogger(MiXCRWorkflowClient.class.getName());
    private String queue;
    private Map<String, SqwFile> tempFiles;

    // meta-types
    private final static String TXT_METATYPE = "text/plain";
    private static final String FASTQ_GZIP_MIMETYPE = "chemical/seq-na-fastq-gzip";

    private void init() {
        try {
            //dir
            dataDir = "data";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            // input samples 
            read1Fastq = getProperty("input_read1_fastq");
            read2Fastq = getProperty("input_read2_fastq");

            //Ext id
            outputFilenamePrefix = getProperty("external_name");

            //Programs
            mixcr = getProperty("mixcr");
            javaHome = getProperty("java_home");

            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            // mixcr
            mixcrMem = Integer.parseInt(getProperty("mixcr_mem"));
            exports = "export LD_LIBRARY_PATH=" + this.javaHome + "/lib" + ":$LD_LIBRARY_PATH" + ";"
                    + "export LD_LIBRARY_PATH=" + this.javaHome + "/jre/lib/amd64/server" + ":$LD_LIBRARY_PATH" + ";"
                    + "export PATH=" + this.javaHome + "/bin" + ":$PATH" + ";";

        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void setupDirectory() {
        init();
        this.addDirectory(dataDir);
        this.addDirectory(tmpDir);
        if (!dataDir.endsWith("/")) {
            dataDir += "/";
        }
        if (!tmpDir.endsWith("/")) {
            tmpDir += "/";
        }
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        SqwFile file0 = this.createFile("read1");
        file0.setSourcePath(read1Fastq);
        file0.setType(FASTQ_GZIP_MIMETYPE);
        file0.setIsInput(true);
        SqwFile file1 = this.createFile("read2");
        file1.setSourcePath(read2Fastq);
        file1.setType(FASTQ_GZIP_MIMETYPE);
        file1.setIsInput(true);
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {

        /**
         * Steps for MiXCR
         */
        // provision files (2) -- clones.clns; clones.det.txt; 
        Job parentJob = null;
        this.alignVdjcaFile = this.dataDir + this.outputFilenamePrefix + ".alignments.vdjca";
        this.alignRescued1File = this.dataDir + this.outputFilenamePrefix + ".alignments_rescued_1.vdjca";
        this.alignRescued2File = this.dataDir + this.outputFilenamePrefix + ".alignments_rescued_2.vdjca";
        this.alignRescued2ExtendFile = this.dataDir + this.outputFilenamePrefix + ".alignments_rescued_2_extended.vdjca";
        this.cloneClnsFile = this.dataDir + this.outputFilenamePrefix + ".clones.clns";
        this.cloneDetTxtFile = this.dataDir + this.outputFilenamePrefix + ".clones.det.txt";

        Job vdjcgenes = alignVdjcGenes();

        Job assembly1 = contigAssembly1();
        assembly1.addParent(vdjcgenes);

        Job assembly2 = contigAssembly2();
        assembly2.addParent(assembly1);

        Job extend = extendAlignment();
        extend.addParent(assembly2);

        Job assembleClones = assembleClonoTypes();
        assembleClones.addParent(extend);

        Job exportClones = exportClonoTypes();
        exportClones.addParent(assembleClones);

        // Provision clones.clns, clones.det.txt{{}}
        //clones.clns is binary output file and clones.det.txt is human readable tab-delimited table text file.
        String clonesFile = this.dataDir + this.outputFilenamePrefix + "clones.clns";
        SqwFile clnsFile = createOutputFile(this.cloneClnsFile, TXT_METATYPE, this.manualOutput);
        clnsFile.getAnnotations().put("MiXCR_clones_clns", "MiXCR");
        exportClones.addFile(clnsFile);

        String cloneTableFile = this.dataDir + this.outputFilenamePrefix + "clones.det.txt";
        SqwFile txtFile = createOutputFile(this.cloneDetTxtFile, TXT_METATYPE, this.manualOutput);
        txtFile.getAnnotations().put("MiXCR_clones_det_txt", "MiXCR");
        exportClones.addFile(txtFile);

    }

    private Job alignVdjcGenes() {
        Job vdjcgenes = getWorkflow().createBashJob("vdjcgenes");
        Command cmd = vdjcgenes.getCommand();
        cmd.addArgument(this.exports);
        cmd.addArgument(this.mixcr);
        cmd.addArgument("align -p rna-seq -s hsa -OallowPartialAlignments=true");
        cmd.addArgument(getFiles().get("read1").getProvisionedPath());
        cmd.addArgument(getFiles().get("read2").getProvisionedPath());
        cmd.addArgument(this.alignVdjcaFile);
        vdjcgenes.setMaxMemory(Integer.toString(mixcrMem * 1024));
        vdjcgenes.setQueue(queue);
        return vdjcgenes;

    }

    private Job contigAssembly1() {
        Job assembly1 = getWorkflow().createBashJob("assembly1");
        Command cmd = assembly1.getCommand();
        cmd.addArgument(this.exports);
        cmd.addArgument(this.mixcr);
        cmd.addArgument("assemblePartial");
        cmd.addArgument(this.alignVdjcaFile);
        cmd.addArgument(this.alignRescued1File);
        assembly1.setMaxMemory(Integer.toString(mixcrMem * 1024));
        assembly1.setQueue(queue);
        return assembly1;
    }

    private Job contigAssembly2() {
        Job assembly2 = getWorkflow().createBashJob("assembly2");
        Command cmd = assembly2.getCommand();
        cmd.addArgument(this.exports);
        cmd.addArgument(this.mixcr);
        cmd.addArgument("assemblePartial");
        cmd.addArgument(this.alignRescued1File);
        cmd.addArgument(this.alignRescued2File);
        assembly2.setMaxMemory(Integer.toString(mixcrMem * 1024));
        assembly2.setQueue(queue);
        return assembly2;
    }

    private Job extendAlignment() {
        Job extend = getWorkflow().createBashJob("extend");
        Command cmd = extend.getCommand();
        cmd.addArgument(this.exports);
        cmd.addArgument(this.mixcr);
        cmd.addArgument("extendAlignments");
        cmd.addArgument(this.alignRescued2File);
        cmd.addArgument(this.alignRescued2ExtendFile);
        extend.setMaxMemory(Integer.toString(mixcrMem * 1024));
        extend.setQueue(queue);
        return extend;
    }

    private Job assembleClonoTypes() {
        Job assembleClones = getWorkflow().createBashJob("assembleClones");
        Command cmd = assembleClones.getCommand();
        cmd.addArgument(this.exports);
        cmd.addArgument(this.mixcr);
        cmd.addArgument("assemble");
        cmd.addArgument(this.alignRescued2ExtendFile);
        cmd.addArgument(this.cloneClnsFile);
        assembleClones.setMaxMemory(Integer.toString(mixcrMem * 1024));
        assembleClones.setQueue(queue);
        return assembleClones;
    }

    private Job exportClonoTypes() {
        Job exportClones = getWorkflow().createBashJob("exportClones");
        Command cmd = exportClones.getCommand();
        cmd.addArgument(this.exports);
        cmd.addArgument(this.mixcr);
        cmd.addArgument("exportClones");
        cmd.addArgument(this.cloneClnsFile);
        cmd.addArgument(this.cloneDetTxtFile);
        exportClones.setMaxMemory(Integer.toString(mixcrMem * 1024));
        exportClones.setQueue(queue);
        return exportClones;
    }

}
