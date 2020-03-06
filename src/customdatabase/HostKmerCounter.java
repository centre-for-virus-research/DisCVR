package customdatabase;

import utilities.PermutationFiles;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

public class HostKmerCounter {

    private String hostFastaFileName;
    private int kmerSize;
    private String outputDirectory;
    private String kAnalyzeOutputDirectory;
    private String kAnalyzeDirectory;
    private int permutationSize;
    private String kAnalyzeOutputFileName;
    private String outputFilePrefix;




    public HostKmerCounter(String hostFastaFileName, int kmerSize, String outputDirectory,
                           String kAnalyzeOutputDirectory, String kAnalyzeDirectory, int permutationSize)
    {
        //Directly passed instance variables.
        this.hostFastaFileName = hostFastaFileName;
        this.kmerSize = kmerSize;
        this.outputDirectory = outputDirectory;
        this.kAnalyzeOutputDirectory = kAnalyzeOutputDirectory;
        this.kAnalyzeDirectory = kAnalyzeDirectory;
        this.permutationSize = permutationSize;
        // Derived instance variables.
        // Exact filename to provide to KAnalyze as a command-line argument.
        this.kAnalyzeOutputFileName = kAnalyzeOutputDirectory + "hostKmers_" + kmerSize;
        // Universal prefix for final output files containing kmer subsets (e.g. temp/hKmers_AAAAA)
        this.outputFilePrefix = outputDirectory + "hKmers_";

    }

    // Runs kAnalyze on host fasta file. kAnalyze output is saved to disk for further processing.
    public void countHostKmers()
    {
        Runtime rt = Runtime.getRuntime();
        Process proc;
        int interVal = -1;

        // Create shell command for KAnalyze. Created as String array for compatibility with Unix shells.
        String kAnalyzeCommand[] = {"java", "-jar", "Xmx3g", kAnalyzeDirectory+File.separator+"kanalyze.jar","count",
                "-k", kmerSize+"", "-o", kAnalyzeOutputFileName, "-f", "fasta", hostFastaFileName, "-rcanonical"};
        // Run kAnalyze
        try {
            proc = rt.exec(kAnalyzeCommand);
            // Wait for the command to complete.
            interVal = proc.waitFor();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        if (interVal != 0)
            System.out.println("Host K-mers counting encountered some errors: "+interVal);

    }

    private String[] createPermutationFiles(){
        //Create and populate array of permutations for given permutation size.
        char set[] = {'A', 'C', 'G', 'T'};
        PermutationFiles PF = new PermutationFiles(permutationSize,set.length);
        String [] permutations = PF.getPermsArray();
        PF.printAllKLength(set,permutationSize,permutations);

        // Write permutations to files with universal prefix (e.g. temp/hKmers_)
        for(int i=0; i<permutations.length;i++){
            String permsFile = outputFilePrefix+(permutations[i]);
            try {
                PrintWriter pw =  new PrintWriter(new BufferedWriter(new FileWriter(permsFile)));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return permutations;
    }

    private void writeHostKmersToPermutationFiles(String firstPerm, String lastPerm){

        int permSize = firstPerm.length();
        long recordSize=0;
        List<String> records = new ArrayList<String>();
        String perm = firstPerm;

        try (BufferedReader br = new BufferedReader(new FileReader(kAnalyzeOutputFileName))){
            for (String line = null; (line = br.readLine()) != null;) {
                line += "\n";
                byte[] buffer = line.getBytes();

                String kmerPerm = line.substring(0, permSize);

                //if the k-mer's perm is the same as perm
                if(kmerPerm.equals(perm)){
                    records.add(line);
                    recordSize += buffer.length;

                }
                else{ //the k-mer's perm is different
                    //write the buffer to a perm file

                    File permFile = new File(outputDirectory+"/hKmers_"+perm);

                    write(permFile, records,recordSize);
                    perm = kmerPerm;
                    recordSize =0;

                    //new list of records and add current k-mer to it
                    records = new ArrayList<String>();
                    records.add(line);
                    recordSize += buffer.length;
                }
            }

            //write the records for the last permFile
            File permFile = new File(outputDirectory+"/hKmers_"+perm);

            write(permFile, records,recordSize);
        } catch (FileNotFoundException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        } catch (IOException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        }

    }

    private void write(File file, List<String> records, long recordSize){

        try {
            final RandomAccessFile raf = new RandomAccessFile(file,"rw");
            raf.seek(raf.length());
            final FileChannel fc = raf.getChannel();
            final MappedByteBuffer mbf = fc.map(FileChannel.MapMode.READ_WRITE,fc.position(), recordSize);
            fc.close();
            for(int i=0;i<records.size();i++){
                final byte [] recordBytes = records.get(i).getBytes(Charset.forName("ISO-8859-1"));
                mbf.put(recordBytes);
            }
            raf.close();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public static void main(String[] args){
        if(args.length < 3)
        {
            printArgumentErrorMessage(args.length);
            System.exit(0);
        }
        else
        {
            printArguments(args);
            String hostFastaFilename = args[0];
            int kmerSize = Integer.parseInt(args[1]);
            String outputDirectoryName = args[2];

            String workingDirectory = System.getProperty("user.dir");
            String kAnalyzeDirectory = workingDirectory + File.separator + "lib";
            int permutationSize = 5;
            outputDirectoryName = workingDirectory + File.separator + outputDirectoryName;
            String kAnalyzeOutputDirectoryName = outputDirectoryName + File.separator + "hostKmers"+kmerSize;

            // Make output directories if they do not exist.
            File outputDirectory = new File(outputDirectoryName);
            if(! outputDirectory.exists()) outputDirectory.mkdir();
            File kAnalyzeOutputDirectory = new File(kAnalyzeOutputDirectoryName);
            if(! kAnalyzeOutputDirectory.exists()) kAnalyzeOutputDirectory.mkdir();

            HostKmerCounter hostKmerCounter= new HostKmerCounter(hostFastaFilename, kmerSize, outputDirectoryName,
                    kAnalyzeOutputDirectoryName, kAnalyzeDirectory, permutationSize);
            hostKmerCounter.countHostKmers();

            String[] hostPermsFiles = hostKmerCounter.createPermutationFiles();
            hostKmerCounter.writeHostKmersToPermutationFiles(hostPermsFiles[0], hostPermsFiles[hostPermsFiles.length-1]);

        }
    }

    private static void printArgumentErrorMessage(int numArguments)
    {
        System.out.println("=================================================================");
        System.out.println();
        System.out.println("Insufficient Number of arguments: "+ numArguments);
        System.out.println("");
        System.out.println("Reguired parameters:");
        System.out.println("[1]: Host genome file");
        System.out.println("[2]: k-mer size");
        System.out.println("[3]: Output directory to store host k-mers in");
        System.out.println();
        System.out.println("=================================================================");
    }

    private static void printArguments(String[] args)
    {
        System.out.println("=================================================================");
        System.out.println();
        System.out.println("Host genome file: "+args[0]);
        System.out.println("k-mer size: "+args[1]);
        System.out.println("Output directory to store host k-mers in: "+args[2]);


        System.out.println("==============================================");
    }



}
