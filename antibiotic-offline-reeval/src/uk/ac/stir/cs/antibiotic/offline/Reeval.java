package uk.ac.stir.cs.antibiotic.offline;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import uk.ac.stir.cs.antibiotic.AntibioticModel;
import uk.ac.stir.cs.antibiotic.AntibioticObjective;
import uk.ac.stir.cs.antibiotic.AntibioticProblem;

public class Reeval {

    public static void main(String[] args) throws IOException, NoSuchAlgorithmException {
        
        Path root = new File(extract(args, "root")).toPath();
        int pid = Integer.parseInt(extract(args, "pid"));
        int mod = Integer.parseInt(extract(args, "mod"));
        int runs = Integer.parseInt(extract(args, "runs"));
        
        if (pid < 0) {
            throw new IllegalArgumentException("pid = " + pid + ", expected: >= 0");
        }
        if (mod < 1) {
            throw new IllegalArgumentException("mod = " + mod + ", expected: >= 1");
        }
        if (pid >= mod) {
            throw new IllegalArgumentException("pid = " + pid + ", expected: < mod");   
        }
        if (runs < 1) {
            throw new IllegalArgumentException("runs = " + runs + ", expected: >= 1");
        }
        
        System.out.println("root = " + root);
        System.out.println("pid = " + pid);
        System.out.println("mod = " + mod);
        System.out.println("runs = " + runs);
        
        //Path root = new File("D:\\Box Sync\\Research\\Antibiotic Resistance Function\\2018 Result\\Full").toPath();
        
        List<Path> vars = Files.walk(root).filter(p-> p.toString().endsWith("VAR.tsv")).collect(Collectors.toList());
                
        AntibioticModel model = AntibioticModel.fixedSampleSize(runs, 700, 100);
        AntibioticProblem problem = 
                AntibioticProblem.builder(model,
                                          AntibioticObjective.uncuredProportion())
                                 .setMaximumIndividualDosage(60)
                                 .build();
        
        for (Path var : vars) {
            File parent = var.getParent().toFile();
            Path fun = new File(parent, "FUN.tsv").toPath();
            if (fun.toFile().exists()) {
                process(parent, pid, mod, problem);
            }
        }
        
    }

    private static void process(File directory, int pid, int modulus, AntibioticProblem problem) throws NoSuchAlgorithmException, IOException {
        File var = new File(directory, "VAR.tsv");
        //File fun = new File(directory, "FUN.tsv");
        File reevaled = new File(directory, "FUN-reevaled.tsv");
        int d = digest(directory, modulus);
        if (d == pid) {
            System.out.println("processing " + directory);
            process(var, /*fun, */reevaled, problem);
        } else {
            System.out.println("skipping " + directory);
        }
    }

    private static int digest(File directory, int modulus) throws NoSuchAlgorithmException {
        MessageDigest md = MessageDigest.getInstance("SHA-256");
        md.update(directory.toString().getBytes());
        byte[] digest = md.digest();
        int last = Byte.toUnsignedInt(digest[digest.length-1]);
        return last % modulus;
    }

    private static void process(File var, /*File fun,*/ File reevaled, AntibioticProblem problem) throws IOException {
        List<String> varLines = lines(var);
//        List<String> funLines = lines(fun);
//        if (varLines.size() != funLines.size()) {
//            throw new RuntimeException("Mismatched line counts for " + var.getParent());
//        }
        Map<String, Double> seen = new HashMap<>();
        int counter = 0;
        try (PrintWriter out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new FileOutputStream(reevaled))))) {
            for (int i = 0; i < varLines.size(); i++) {
                System.out.print("line " + counter + "\t");
                String varLine = varLines.get(i);
//                String funLine = funLines.get(i);
                int[] x = parseCandidate(varLine);
//                double[] f = parsePair(funLine);
                double fitness;
                if (seen.containsKey(varLine)) {
                    fitness = seen.get(varLine);
                    System.out.println("skipped");
                } else {
                    fitness = problem.evaluateFirstObjective(x);
                    seen.put(varLine, fitness);
                    System.out.println("processed");
                }
                out.println(fitness);
                out.flush();
                counter++;
            }
            out.flush();
        }
    }

    private static List<String> lines(File file) throws IOException {
        List<String> lines = new ArrayList<>(100);
        String line;
        try (BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(file)))) {
            while ((line = in.readLine()) != null) {
                if (line.trim().length() > 0) {
                    lines.add(line.trim());
                }
            }    
        }
        return lines;
    }

    private static int[] parseCandidate(String line) {
        String[] tokens = line.split(" ");
        if (tokens.length > 10) {
            throw new RuntimeException("line had wrong length : " + line);
        }
        int[] rv = new int[10];
        for (int i = 0; i < tokens.length; i++) {
            rv[i] = Integer.parseInt(tokens[i]);
        }
        return rv;
    }

//    private static double[] parsePair(String line) {
//        String[] tokens = line.split(" ");
//        if (tokens.length != 2) {
//            throw new RuntimeException("line had wrong length : " + line);
//        }
//        double[] rv = new double[2];
//        for (int i = 0; i < 2; i++) {
//            rv[i] = Double.parseDouble(tokens[i]);
//        }
//        return rv;
//    }

    private static String extract(String[] args, String name) {
        int start = indexOf(name, args);
        return args[start+1];
    }
    
//    private static String extractQuoted(String[] args, String name) {
//        int start = indexOf(name, args);
//        if (!args[start+1].startsWith("\"")) {
//            System.out.println(args[start+1]);
//            throw new RuntimeException("Missing \" for " + name);
//        }
//        StringBuilder rv = new StringBuilder();
//        for (int i = start + 1; i < args.length; i++) {
//            if (args[i].startsWith("--")) {
//                break;
//            }
//            if (args[i].startsWith("\"")) {
//                rv.append(args[i].substring(1));
//            } else {
//                rv.append(args[i]);
//            }
//            if (args[i].endsWith("\"")) {
//                break;
//            }
//        }
//        return rv.substring(0, rv.length()-1);
//    }

    private static int indexOf(String name, String[] args) {
        for (int i = 0; i < args.length; i++) {
            if (args[i].equalsIgnoreCase("--" + name)) {
                return i;
            }
        }
        throw new RuntimeException("--" + name + " not found in args");
    }
    
}
