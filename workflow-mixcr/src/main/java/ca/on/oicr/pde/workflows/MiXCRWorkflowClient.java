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
 * setupFiles 3. setupWorkflow 4. setupEnvironment 5. buildWorkflow
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
    private String alignvdjcaFile;
    private String alignrescued1File;
    private String alignrescued2File;
    private String alignrescued2extendFile;
    private String cloneClnsFile;
    private String cloneDetTxtFile;
    
    
    //Tools
    private String Mixcr;
    private String java;


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
   // private final static String TAR_GZ_METATYPE = "application/tar-gzip";
   // private static final String FASTQ_GZIP_MIMETYPE = "chemical/seq-na-fastq-gzip";

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
            Mixcr = getProperty("MIXCR");
            java = getProperty("JAVA");
            
            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            // mixcr
            mixcrMem = Integer.parseInt(getProperty("mixcr_mem"));
   


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
        file0.setType(TXT_METATYPE);
        file0.setIsInput(true);
        SqwFile file1 = this.createFile("read2");
        file1.setSourcePath(read2Fastq);
        file1.setType(TXT_METATYPE);
        file1.setIsInput(true);
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {

        /**
         * Steps for sequenza: 1. Check if "bam" file exists; true 2. Check if
         * "bai" file exists; true: go to step 4 3. Check if normal Pb_R sample
         * exists; true: go to step 4; else abort 3. If false: samtools index
         * "bam" file 4. Run job sequenza-utils 5. If outputFile ends with
         * "bin50.gz"; go to step 6; else go to step 4 6. Run job sequenzaR 7.
         * Iterate through the files/folders in outDir: 8. If fileName1 ==
         * "pandc.txt" and fileName2 ends with "Total_CN.seg"; create a folder
         * called "copynumber" 9. If fileType == "folder"; create a folder
         * called "model-fit"; move folders to "model-fit" 10. If fileType ==
         * "file" && fileName != outputFile; move file to "model-fit" 11. Delete
         * outputFile (rm outputFile) 12. zip "model-fit" 13. outputFile =
         * fileName2 14. OutputDir contains the following: fileName1,
         * outputFile, model-fit.zip
         */
        // workflow : read inputs tumor and normal bam; run sequenza-utils; write the output to temp directory; 
        // run sequenzaR; handle output; provision files (3) -- model-fit.zip; text/plain; text/plain
        
        Job parentJob = null;
        this.alignvdjcaFile =this.dataDir+  this.outputFilenamePrefix + ".alignments.vdjca";
        this.alignrescued1File = this.dataDir + this.outputFilenamePrefix + ".alignments_rescued_1.vdjca";
        this.alignrescued2File = this.dataDir + this.outputFilenamePrefix + ".alignments_rescued_2.vdjca";
        this.alignrescued2extendFile = this.dataDir + this.outputFilenamePrefix + ".alignments_rescued_2_extended.vdjca";
        this.cloneClnsFile = this.dataDir + this.outputFilenamePrefix + ".clones.clns";
        this.cloneDetTxtFile = this.dataDir + this.outputFilenamePrefix + ".clones.det.txt";
        

        Job VDJCgenes = alignVDJCgenes();
        parentJob = VDJCgenes;

        Job Assembly1= contigAssembly1();
        Assembly1.addParent(parentJob);
        parentJob = Assembly1;

      
        Job Assembly2 = contigAssembly2();
        Assembly2.addParent(parentJob);
        parentJob = Assembly2;

        Job Extend = extendAlignment();
        Extend .addParent(parentJob);
        parentJob = Extend ;

        Job AssembleClones = assembleClonotypes();
        AssembleClones.addParent(parentJob);
        parentJob = AssembleClones;

        Job ExportClones = exportClonotypes();
        ExportClones.addParent(parentJob);

        // Provision clones.clns, clones.det.txt{{}}
        String clonesFile = this.tmpDir + "clones.clns";
        SqwFile clnsFile = createOutputFile(clonesFile, TXT_METATYPE, this.manualOutput);
        clnsFile.getAnnotations().put("MiXCR_clones_clns", "MiXCR");
        ExportClones.addFile(clnsFile);
        
        String cloneTableFile = this.tmpDir + "clones.det.txt";
        SqwFile txtFile = createOutputFile(cloneTableFile, TXT_METATYPE, this.manualOutput);
        txtFile.getAnnotations().put("MiXCR_clones_det_txt", "MiXCR");
        ExportClones.addFile(txtFile);  
    
    }
    
    
    
    private Job alignVDJCgenes() {
        Job VDJCgenes = getWorkflow().createBashJob("VDJCgenes");
        Command cmd = VDJCgenes.getCommand();
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/lib" + ";");
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/jre/lib/amd64/server" + ";");
        cmd.addArgument("export PATH=" + this.java + "/bin" + ";");
        cmd.addArgument(this.Mixcr);
        cmd.addArgument("align -p rna-seq -s hsa -OallowPartialAlignments=true");
        cmd.addArgument(getFiles().get("read1").getProvisionedPath());
        cmd.addArgument(getFiles().get("read2").getProvisionedPath());
        cmd.addArgument(this.alignvdjcaFile);
        VDJCgenes.setMaxMemory(Integer.toString(mixcrMem * 1024));
        VDJCgenes.setQueue(getOptionalProperty("queue", ""));
        return VDJCgenes;
        
    }
    
      private Job contigAssembly1() {
        Job Assembly1 = getWorkflow().createBashJob("Assembly1");
        Command cmd = Assembly1.getCommand();
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/lib" + ";");
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/jre/lib/amd64/server" + ";");
        cmd.addArgument("export PATH=" + this.java + "/bin" + ";");
        cmd.addArgument(this.Mixcr);
        cmd.addArgument("assemblePartial");
        cmd.addArgument(this.alignvdjcaFile);
        cmd.addArgument(this.alignrescued1File);
        Assembly1.setMaxMemory(Integer.toString(mixcrMem * 1024));
        Assembly1.setQueue(getOptionalProperty("queue", ""));
        return Assembly1;
    }
    
      
      private Job contigAssembly2() {
        Job Assembly2 = getWorkflow().createBashJob("Assembly2");
        Command cmd = Assembly2.getCommand();
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/lib" + ";");
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/jre/lib/amd64/server" + ";");
        cmd.addArgument("export PATH=" + this.java + "/bin" + ";");
        cmd.addArgument(this.Mixcr);
        cmd.addArgument("assemblePartial");
        cmd.addArgument(this.alignrescued1File);
        cmd.addArgument(this.alignrescued2File);
        Assembly2.setMaxMemory(Integer.toString(mixcrMem * 1024));
        Assembly2.setQueue(getOptionalProperty("queue", ""));
        return Assembly2;
    }
    
      
      private Job extendAlignment() {
        Job Extend = getWorkflow().createBashJob("Extend");
        Command cmd = Extend.getCommand();
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/lib" + ";");
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/jre/lib/amd64/server" + ";");
        cmd.addArgument("export PATH=" + this.java + "/bin" + ";");
        cmd.addArgument(this.Mixcr);
        cmd.addArgument("extendAlignments");
        cmd.addArgument(this.alignrescued2File);
        cmd.addArgument(this.alignrescued2extendFile);
        Extend.setMaxMemory(Integer.toString(mixcrMem * 1024));
        Extend.setQueue(getOptionalProperty("queue", ""));
        return Extend;
    }   
      
      
     private Job assembleClonotypes() {
        Job AssembleClones = getWorkflow().createBashJob("AssembleClones");
        Command cmd = AssembleClones.getCommand();
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/lib" + ";");
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/jre/lib/amd64/server" + ";");
        cmd.addArgument("export PATH=" + this.java + "/bin" + ";");
        cmd.addArgument(this.Mixcr);
        cmd.addArgument("assemble");
        cmd.addArgument(this.alignrescued2extendFile);
        cmd.addArgument(this.cloneClnsFile);
        AssembleClones.setMaxMemory(Integer.toString(mixcrMem * 1024));
       AssembleClones.setQueue(getOptionalProperty("queue", ""));
        return AssembleClones;
    }   
     
     private Job exportClonotypes() {
        Job ExportClones = getWorkflow().createBashJob("ExportClones");
        Command cmd = ExportClones.getCommand();
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/lib" + ";");
        cmd.addArgument("export LD_LIBRARY_PATH=" + this.java + "/jre/lib/amd64/server" + ";");
        cmd.addArgument("export PATH=" + this.java + "/bin" + ";");
        cmd.addArgument(this.Mixcr);
        cmd.addArgument("exportClones");
        cmd.addArgument(this.cloneClnsFile);
        cmd.addArgument(this.cloneDetTxtFile);
        ExportClones.setMaxMemory(Integer.toString(mixcrMem * 1024));
        ExportClones.setQueue(getOptionalProperty("queue", ""));
        return ExportClones;
    }   
       
       
      
}
