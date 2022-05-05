/**
 * Created by jayna on 17/05/2018.
 */

import cern.colt.list.DoubleArrayList;
import cern.jet.random.Poisson;
import cern.jet.random.Uniform;
import mpi.MPI;
import mpi.MPIException;
import mpi.Status;
import org.apache.commons.lang3.ArrayUtils;

import java.io.*;
import java.util.*;
import org.nustaq.serialization.*;
import org.nustaq.serialization.simpleapi.DefaultCoder;


public class patchSimMPI {


    params params;
    infectionHistory patch_history;

    double tau; // time step
    int genotype;
    Map<Integer, List<Integer>> genotype_curr_per_patch;
    List<Double> potentialSamplingTimepoints;
    List<Integer> emptyPatchList;
    List<Integer> occupyPatchList;
    int total_infected;
    int total_genotypes;
    List<Integer> X_curr_i;
    List<Integer> Y_curr_i;
    List<Integer> curr_in_body;
    int [] threshold_i;
    List<Double> timeThreshold_curr_i;
    List<Double>  timeRefactory_curr_i;
    static FSTConfiguration conf;

    public patchSimMPI() {

        tau = 0.25;
        patch_history = new infectionHistory();
        genotype = 1;
        genotype_curr_per_patch = new HashMap<>();
        emptyPatchList = new ArrayList<>();
        occupyPatchList = new ArrayList<>();
        total_genotypes = 0;
        total_infected = 0;
        X_curr_i = new ArrayList<>();
        Y_curr_i = new ArrayList<>();
        curr_in_body = new ArrayList<>();
        conf = FSTConfiguration.createDefaultConfiguration();


    }

    public void runSim(params params) throws MPIException{

        /* MPI variables */

        int rank, size, source, tag;

        /* Other variables */

        int istart = 0, istop = 0;

        size = MPI.COMM_WORLD.getSize();
        rank = MPI.COMM_WORLD.getRank();

        String outputfile = "patchSim_extinction_"+params.nu+"_col_"+params.c+
                "_beta_"+params.beta+"_Npatches_"+params.Npatches+
                "_patchSize_"+(params.S+params.I+".txt");


        String summaryFile = "summary_patchSim_extinction_"+params.nu+"_col_"+params.c+
                "_beta_"+params.beta+"_Npatches_"+params.Npatches+
                "_patchSize_"+(params.S+params.I+".txt");

        FileWriter writer1 = null;
        FileWriter writer2 = null;

        try {
            writer1 = new FileWriter(new File(outputfile));
            writer1.write("Time\tPatch_no\tTotal_infected\tUnique_genotypes\tGenealogical_diversity\tTMRCA\tGenetic_divergence\tSim\n");

            writer2 = new FileWriter(new File(summaryFile));
            writer2.write("Time\tTotal_infected\tUnique_genotypes\tOccupied_patches\tGenealogical_diversity\tTMRCA\tGenetic_divergence\tSim\n");
        } catch (IOException e) {
            e.printStackTrace();
        }


        double t_curr = tau;
        double t_max = params.runTime;
        int Y_curr = params.I;
        int X_curr = params.S;
        DoubleArrayList withinPatchRates = new DoubleArrayList();
        DoubleArrayList betweenPatchRates = new DoubleArrayList();


        double[] beta_i = new double[params.Npatches];
        int[] daysRefactory_i = new int[params.Npatches];


        //initializing patches with number of infected and susceptible cells


//        Arrays.fill(timeThreshold_curr_i, 0.0);
//        Arrays.fill(timeRefactory_curr_i, 0);
        Arrays.fill(beta_i, params.beta);
        Arrays.fill(daysRefactory_i, 5);
        Y_curr_i = new ArrayList<>(Collections.nCopies(params.Npatches, 0));
        X_curr_i = new ArrayList<>(Collections.nCopies(params.Npatches, params.S+params.I));
        timeThreshold_curr_i = new ArrayList<>(Collections.nCopies(params.Npatches, 0.0));
        timeRefactory_curr_i = new ArrayList<>(Collections.nCopies(params.Npatches, 0.0));

        threshold_i = new int[params.Npatches];

        double colonize = 0.0;
        double extinct = 0.0;
        double death = 0.1;


        Y_curr_i.set(0, Y_curr);
        X_curr_i.set(0, X_curr_i.get(0)-Y_curr);

        initialiseHistory(patch_history, Y_curr_i, t_max, 1);

        for (int p = 0; p < params.Npatches; p++) {

            List<Integer> genotype_curr_i = new ArrayList<>();

            for(Integer geno: patch_history.genotype) {

                int genotype_index = patch_history.genotype.indexOf(geno);

                if(patch_history.patch.get(genotype_index).equals(p)) {

                    int Y_temp = patch_history.prevalence.get(genotype_index).get(0);

                    int s = 0;
                    while(s < Y_temp) {
                        genotype_curr_i.add(geno);
                        s++;
                    }
                }
            }

            genotype_curr_per_patch.put(p, genotype_curr_i);

        }

        //controls the extinction rate...


        for(int i = 0; i < threshold_i.length; i++) {

            int N = params.S+params.I;
            threshold_i[i] = Uniform.staticNextIntFromTo((int)Math.ceil(0.5*N),(int)Math.ceil(0.9*N));

            if(genotype_curr_per_patch.get(i).size() == 0) {
                emptyPatchList.add(i);
            }
            else{
                occupyPatchList.add(i);
            }
        }


        byte [] patch_bytes = null;

        //Send current number to the right and receive from the left
        int left = rank -1;
        int right = rank + 1;
        if ((rank+1) == size){
            right = 0;
        }
        if(rank == 0){
            left = size - 1;
        }

        int[] test = new int[1];
        test[0] = 1;


        if(size == 2) {

            if(rank == 1) {
                System.out.println("x, "+size);
                istart = (params.Npatches / size) * (rank - 1);
                istop = istart + (params.Npatches);


            }
        }
        else if(size > 2) {

            //System.out.println("y, "+size+", "+rank);

            istart = ((params.Npatches / (size-1)) * (rank-1));
            istop = istart + (params.Npatches/(size-1));

        }

        int sum = 0;

        DefaultCoder coder = null;

        if(rank == 0) {
            System.out.println(params.Npatches+", "+params.c+", "+params.nu+", "+params.runTime);
        }

        for (int t = 1; t < (t_max / tau); t++) {
            t_curr += tau;

            rank = MPI.COMM_WORLD.getRank();
            Poisson poisson;
            int j;
            if(rank == 0) {
                curr_in_body.clear();

                String sim = "ext_" + extinct + "_col_" + colonize + "_Npatch_" + params.Npatches + "_PatchSize_" + (params.S + params.I);


                updatePatchStatistics(writer1, t_curr, sim);

                //determine what is present currently across all patches
                update_curr_in_body(curr_in_body, t);


                //patch_bytes = serializeData();

                if (t_curr % 10 == 0) {

                    List<Double> globalDiversity = updateDiversity(curr_in_body, true);

                    //double geneticDiversity = globalDiversity.get(2);
                    double divergence = globalDiversity.get(2);


                    System.out.println("t_curr \t" + t_curr + "\ttotal_infected: \t" + Math.log10(total_infected) + "\ttotal_genotypes: \t" + total_genotypes +
                            "\tpatch_occupied: \t" + occupyPatchList.size() + "\t" + divergence);


                    writeOutput(writer2, t_curr, total_infected, total_genotypes, occupyPatchList.size(), globalDiversity, sim);


                }

                colonize = params.c;
                extinct = params.nu;

                if (t_curr < 40) {

                    extinct = 0;
                    colonize = 0.1;

                }

                betweenPatchRates.add(extinct * occupyPatchList.size());
                betweenPatchRates.add(colonize * emptyPatchList.size());

                for (int event = 0; event < betweenPatchRates.size(); event++) {

                    poisson = new Poisson(tau * betweenPatchRates.get(event), params.randomGenerator);
                    int num = poisson.nextInt();

                    if (event == 0) {

                        j = 0;

                        //extinction

                        int min = Math.min(num, occupyPatchList.size());

                        while (j < min) {

                            // patch extinction

                            int r = Uniform.staticNextIntFromTo(0, occupyPatchList.size() - 1);
                            int patch_index = occupyPatchList.get(r);


                            //if (Y_curr_i.get(patch_index) >= threshold_i[patch_index]) {


                            Y_curr_i.set(patch_index, 0);
                            X_curr_i.set(patch_index, params.S + params.I);

                            //randomly choose a genotype to seed/colonize new patch

                            timeThreshold_curr_i.set(patch_index, 0.0);
                            timeRefactory_curr_i.set(patch_index, t_curr);
                            // track infection history (genotype_curr should be set to 0)

                            Set<Integer> uniqueGenotypes = new HashSet<>();
                            uniqueGenotypes.addAll(genotype_curr_per_patch.get(patch_index));

                            for (Integer g : uniqueGenotypes) {

                                int index = patch_history.genotype.indexOf(g);

                                if (patch_history.patch.get(index) == patch_index) {
                                    patch_history.prevalence.get(index).set(t, 0);
                                    patch_history.death.set(index, t_curr);

                                }

                            }

                            //threshold_i[patch_index] = Uniform.staticNextIntFromTo((int) Math.ceil(0.7 * (params.S + params.I)), (int) Math.ceil(0.9 * (params.S + params.I)));

                            genotype_curr_per_patch.get(patch_index).clear();

                            occupyPatchList.remove((Integer) patch_index);
                            emptyPatchList.add(patch_index);

                            //}

                            j++;
                        }

                    }
                    else {

                        // patch colonisation

                        j = 0;

                        int min = Math.min(num, emptyPatchList.size());

                        while (j < min) {

                            int r = Uniform.staticNextIntFromTo(0, emptyPatchList.size() - 1);

                            int patch_index = emptyPatchList.get(r);

                            int probTimeRefactory = daysRefactory_i[patch_index]; // choose probabilistically how many days patch remains in refactory period (make a variable of the model)

                            if ((t_curr - timeRefactory_curr_i.get(patch_index)) > (double) probTimeRefactory) {

                                if (curr_in_body.size() == 0) {

                                    break;
                                }

//                            if (curr_in_body.size() == 0) {
//
//                                while(curr_in_body.size() == 0) {
//
//                                    double delta = Uniform.staticNextDoubleFromTo(0, 2000);
//
//                                    if(delta >= t_curr) {
//
//                                        delta = t_curr;
//                                    }
//                                    double chosen_time_in_the_past = Uniform.staticNextDoubleFromTo(t_curr-delta, t_curr);
//
//
//                                    for (int i = 0; i < patch_history.genotype.size(); i++) {
//
//                                        if (patch_history.birth.get(i) <= chosen_time_in_the_past && patch_history.death.get(i) >= chosen_time_in_the_past) {
//
//                                            curr_in_body.add(patch_history.genotype.get(i));
//
//                                        }
//                                    }
//                                }
//                            }

                                Collections.shuffle(curr_in_body);

                                timeRefactory_curr_i.set(patch_index, 0.0);

                                int newInfections = Uniform.staticNextIntFromTo(1, 3);

                                Y_curr_i.set(patch_index, newInfections);
                                X_curr_i.set(patch_index, (X_curr + Y_curr) - newInfections);

                                r = Uniform.staticNextIntFromTo(0, curr_in_body.size() - 1);
                                Integer colonizing_genotype = curr_in_body.get(r);


                                for (int i = 0; i < newInfections; i++) {

                                    genotype_curr_per_patch.get(patch_index).add(genotype);

                                }
                                patch_history.logBirth(t_curr);
                                patch_history.logDeath(Double.POSITIVE_INFINITY);
                                patch_history.logGenotype(genotype);
                                patch_history.logMutations(0);
                                patch_history.logParent(colonizing_genotype);
                                patch_history.logPatch(patch_index);

                                List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau)]));
                                Collections.fill(new_prevalence, 0);
                                new_prevalence.set(t, newInfections);
                                patch_history.prevalence.add(new_prevalence);
                                genotype++;

                                emptyPatchList.remove((Integer) patch_index);
                                occupyPatchList.add((Integer) (patch_index));

                            }
                            j++;
                        }


                    }

                }
                betweenPatchRates.removeAll(betweenPatchRates);


                if(coder == null) {
                    coder = new DefaultCoder();
                }
                patch_bytes = serializeData(coder);

                int bytes_length = Math.round(patch_bytes.length);
                int sub_length = (int)Math.ceil((double)bytes_length/6.0);


                for(int s=0; s < 5; s++) {
                    byte[] temp = Arrays.copyOfRange(patch_bytes, sub_length*s, sub_length*(s+1));
                    MPI.COMM_WORLD.send(temp, temp.length, MPI.BYTE, right, (s+1));

                }
                byte[] temp = Arrays.copyOfRange(patch_bytes, sub_length*5, bytes_length);
                MPI.COMM_WORLD.send(temp, temp.length, MPI.BYTE, right, (5+1));

            }
            else {

                byte[] x = null;
                for(int s=0; s < 6; s++) {
                    Status status = MPI.COMM_WORLD.probe(left, (s+1));
                    int m_length = status.getCount(MPI.BYTE);
                    byte[] recv = new byte[m_length];
                    MPI.COMM_WORLD.recv(recv, recv.length, MPI.BYTE, left, (s+1));

                    x = ArrayUtils.addAll(x, recv);

                }

                if(coder == null) {
                    coder = new DefaultCoder();
                }
                deserializaData(x, coder);


                for (int p = istart; p < istop; p++) {

                    int X_temp = X_curr_i.get(p);
                    int Y_temp = Y_curr_i.get(p);

                    withinPatchRates.add((beta_i[p] * X_temp * Y_temp) / (X_temp + Y_temp)); // infection rate
                    withinPatchRates.add(death * (Y_temp));

                    //Poisson poisson;
                    int noOfRates = withinPatchRates.size();
                    //int j;


                    //events in patches that affect X and Y individuals
                    for (int event = 0; event < noOfRates; event++) {

                        poisson = new Poisson(tau * withinPatchRates.get(event), params.randomGenerator);
                        int num = poisson.nextInt();

                        //int patch_index = p;

                        List<Integer> genotype_curr = new ArrayList<>(genotype_curr_per_patch.get(p));


                        if (event == 0) {

                            // infection

                            int min = Math.min(X_curr_i.get(p), num);

                            j = 0;


                            while (j < min) {

                                //determine number of mutations
//                                int mutations = Poisson.staticNextInt(params.U); //double check if rate vs mean


//                                if (mutations == 0) {
                                    int r = (int) Math.floor(Uniform.staticNextDouble() * genotype_curr.size());
                                    Integer chosen_genotype = genotype_curr.get(r);
                                    int genotype_index = patch_history.genotype.indexOf(chosen_genotype);

//                                List<Integer> prevalence = new ArrayList<>();
//                                prevalence.addAll(patch_history.prevalence.get(genotype_index));

                                    int prevalence_curr = patch_history.prevalence.get(genotype_index).get(t - 1);
                                    int prevalence_next = patch_history.prevalence.get(genotype_index).get(t);

                                    X_curr_i.set(p, X_curr_i.get(p) - 1);
                                    Y_curr_i.set(p, Y_curr_i.get(p) + 1);

                                    if (Y_curr_i.get(p) >= threshold_i[p]) {

                                        timeThreshold_curr_i.set(p, t_curr);

                                    }

                                    if (prevalence_next > 0) {

                                        prevalence_curr = prevalence_next;
                                    }

                                    genotype_curr.add(chosen_genotype);
                                    patch_history.prevalence.get(genotype_index).set(t, prevalence_curr + 1);
//                                }
//                                else {
//
//
//                                    //choose genotype
//                                    int index = (int) Math.floor(Uniform.staticNextDouble() * genotype_curr.size());
//                                    Integer chosen_genotype = genotype_curr.get(index);
//                                    int parent_index = patch_history.genotype.indexOf(chosen_genotype);
//
//
//                                    patch_history.logParent(chosen_genotype); // do we want to log parent index or genotype?
//                                    patch_history.logGenotype(genotype);
//                                    patch_history.logMutations(mutations);
//                                    patch_history.logBirth(t_curr);
//                                    patch_history.logDeath(Double.POSITIVE_INFINITY);
//                                    patch_history.logPatch(p);
//
////                                List<Integer> prevalence = new ArrayList<>();
////                                prevalence.addAll(patch_history.prevalence.get(parent_index));
//
//                                    //prevalence is bit off need to check this
//                                    int prevalence_curr = patch_history.prevalence.get(parent_index).get(t - 1);
//                                    int prevalence_next = patch_history.prevalence.get(parent_index).get(t);
//
//                                    if (prevalence_next > 0) {
//
//                                        prevalence_curr = prevalence_next;
//                                    }
//
//                                    if ((prevalence_curr - 1) == 0) {
//
//                                        patch_history.death.set(parent_index, t_curr);
//                                    }
//                                    patch_history.prevalence.get(parent_index).set(t, prevalence_curr - 1);
//
//                                    List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[patch_history.prevalence.get(parent_index).size()]));
//                                    Collections.fill(new_prevalence, 0);
//                                    new_prevalence.set(t, 1);
//                                    patch_history.prevalence.add(new_prevalence);
//
//                                    genotype_curr.remove(chosen_genotype);
//                                    genotype_curr.add(genotype);
//                                    genotype++;
//
//                                }
                                j++;
                            }

                        }
                        else {

                            int min = Math.min(Y_curr_i.get(p), num);

                            j = 0;

                            while (j < min) {

                                X_curr_i.set(p, (X_curr_i.get(p) + 1));
                                Y_curr_i.set(p, (Y_curr_i.get(p) - 1));

                                int index = (int) Math.floor(Uniform.staticNextDouble() * genotype_curr.size());
                                Integer chosen_genotype = genotype_curr.get(index);
                                int genotype_index = patch_history.genotype.indexOf(chosen_genotype);

                                int prevalence_curr = patch_history.prevalence.get(genotype_index).get(t);
//                            int prevalence_next = prevalence.get(t);
//
//                            if (prevalence_next > 0) {
//
//                                prevalence_curr = prevalence_next;
//                            }

                                if ((prevalence_curr - 1) == 0) {

                                    patch_history.death.set(genotype_index, t_curr);
                                }

                                patch_history.prevalence.get(genotype_index).set(t, prevalence_curr - 1);
                                genotype_curr.remove(chosen_genotype);
                                j++;
                            }
                        }
                        genotype_curr_per_patch.replace(p, genotype_curr);

                    }

                    withinPatchRates.removeAll(withinPatchRates);


                }

                if(coder == null) {
                    coder = new DefaultCoder();
                }
                patch_bytes = serializeData(coder);
                int bytes_length = Math.round(patch_bytes.length);
                int sub_length = (int)Math.ceil((double)bytes_length/6.0);

                for(int s=0; s < 5; s++) {
                    byte[] send = Arrays.copyOfRange(patch_bytes, sub_length*s, sub_length*(s+1));
                    MPI.COMM_WORLD.send(send, send.length, MPI.BYTE, right, (s+1));

                }
                byte[] temp = Arrays.copyOfRange(patch_bytes, sub_length*5, bytes_length);
                MPI.COMM_WORLD.send(temp, temp.length, MPI.BYTE, right, (5+1));
            }
            if(rank == 0) {

                byte[] x = null;
                for(int s=0; s < 6; s++) {
                    Status status = MPI.COMM_WORLD.probe(left, (s+1));
                    int m_length = status.getCount(MPI.BYTE);
                    byte[] recv = new byte[m_length];
                    MPI.COMM_WORLD.recv(recv, recv.length, MPI.BYTE, left, (s+1));
                    x = ArrayUtils.addAll(x, recv);

                }
                if(coder == null) {
                    coder = new DefaultCoder();
                }
                deserializaData(x, coder);
//                System.out.println(rank+": "+curr_in_body);

            }

        }


        if(rank == 0) {
            try {
                writer1.close();
                writer2.close();

            } catch (IOException e) {
                e.printStackTrace();
            }

        }
    }



    private void initialiseHistory(infectionHistory history, List<Integer> Y_curr_i, double t_max, int n_patches) {

        int i = 0;
        while(i < n_patches) {

            if(Y_curr_i.get(i) > 0) {
                history.logGenotype(genotype);
                history.logMutations(0);
                history.logBirth(Double.NEGATIVE_INFINITY);
                history.logDeath(Double.POSITIVE_INFINITY);
                history.logParent((int) Double.NEGATIVE_INFINITY);
                List<Integer> prevalence_through_time = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau)]));
                Collections.fill(prevalence_through_time, 0);
                prevalence_through_time.set(0, Y_curr_i.get(i));
                history.prevalence.add(prevalence_through_time);
                history.logPatch(i);
                genotype++;
            }
            i++;
        }
    }

    public infectionHistory getInfectionHistory() {

        return this.patch_history;
    }

    private List<Integer> deepCopy(List<Integer> list) {
        List<Integer> copy = new ArrayList<Integer>(list.size());
        for (Integer element : list) {
            Integer elementCopy = element;
            copy.add(elementCopy);
        }
        return copy;
    }

    private void writeOutput(FileWriter writer, double t, int total_infected, int genotypes, int n_patches, List<Double> diversity, String sim) {

        try {
            writer.write(t+"\t"+total_infected+"\t"+genotypes+"\t"+n_patches+"\t"+diversity.get(0)+"\t"+diversity.get(1) + "\t"+diversity.get(2) + "\t"+sim+"\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private int commonAncestor(Integer i1, Integer i2) {

        Integer commonAnc = (int)Double.NEGATIVE_INFINITY;
        Integer lineageA = i1;
        Integer lineageB = i2;


        Set<Integer> ancestry = new HashSet<Integer>();
        while (true) {

            Integer lineage_a = patch_history.getParent(lineageA);
            Integer lineage_b = patch_history.getParent(lineageB);


            if (lineage_a.intValue() != (int) Double.NEGATIVE_INFINITY) {

                lineageA = lineage_a;
                if (!ancestry.add(lineageA)) {
                    commonAnc = lineageA;
                    break;
                }
            }
            if (lineage_b.intValue() != (int) Double.NEGATIVE_INFINITY) {

                lineageB = lineage_b;
                if (!ancestry.add(lineageB)) {
                    commonAnc = lineageB;
                    break;
                }
            }
            if (lineage_a.intValue() == (int) Double.NEGATIVE_INFINITY && lineage_b.intValue() == (int) Double.NEGATIVE_INFINITY) {
                break;
            }

        }


        return commonAnc;
    }

    private double distance(Integer i1, Integer i2) {

        Integer ancestor = commonAncestor(i1,i2);
        if(ancestor.intValue() != (int)Double.NEGATIVE_INFINITY) {

            double ancestorDistance= patch_history.getBirth(ancestor);

            if(Double.isInfinite(ancestorDistance)) {
                ancestorDistance = 0;
            }

            double distA = patch_history.getBirth(i1) - ancestorDistance;
            double distB = patch_history.getBirth(i2) - ancestorDistance;
            return distA + distB;
        }
        else{
            return 0;

        }
    }

    private double geneticDistance(Integer i1, Integer i2) {

        Integer ancestor = commonAncestor(i1,i2);
        //if(ancestor.intValue() != (int)Double.NEGATIVE_INFINITY) {

        double distA = patch_history.getMutationsFromParent(i1, ancestor);
        double distB = patch_history.getMutationsFromParent(i2, ancestor);
        return distA + distB;

    }

    public void updatePatchStatistics(FileWriter writer, double t_curr, String sim) {

        total_infected = 0;
        total_genotypes = 0;

        for(Integer patch: genotype_curr_per_patch.keySet()) {


            total_infected+=genotype_curr_per_patch.get(patch).size();
            total_genotypes+=(new HashSet<>(genotype_curr_per_patch.get(patch))).size();
//            if(genotype_curr_per_patch.get(patch).size()==0) {
//                emptyPatchList.add(patch);
//            }
//            else{
//                occupyPatchList.add(patch);
//            }



            if(t_curr % 10 == 0) {
                List<Double> diversity_results = updateDiversity(genotype_curr_per_patch.get(patch), false);

                double diversity = diversity_results.get(0);
                double tmrca = diversity_results.get(1);
                //double geneticDiversity = diversity_results.get(2);

                try {
                    writer.write(t_curr + "\t" + (patch + 1) + "\t" + genotype_curr_per_patch.get(patch).size() +
                            "\t" + (new HashSet<>(genotype_curr_per_patch.get(patch))).size() +
                            "\t" + diversity + "\t" + tmrca + "\t" + sim + "\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

        }
    }

    public void update_curr_in_body(List<Integer> curr_in_body, int t) {
        for (int i = 0; i < patch_history.prevalence.size(); i++) {

            //List<Integer> prevalence = patch_history.prevalence.get(i);
            int prev_curr = patch_history.prevalence.get(i).get(t - 1);
            patch_history.prevalence.get(i).set(t, prev_curr);


            Integer genotype_i = patch_history.genotype.get(i);

            if (prev_curr > 0) {

//                Integer[] genotype_array = new Integer[patch_history.prevalence.get(i).get(t - 1)];
//                Arrays.fill(genotype_array, genotype_i);

                List<Integer> tmp = Collections.nCopies(prev_curr, genotype_i);
                //curr_in_body.addAll(Arrays.asList(genotype_array));
                curr_in_body.addAll(tmp);
            }

        }

    }

    public List<Double> updateDiversity(List<Integer> genotype_curr, boolean getDivergence) {

        double diversity = 0.0;
        double tmrca1 = 0.0;
        int sampleCount1 = 0;

        //double geneticDiversity = 0.0;
        double divergence = 0.0;

        int limit = genotype_curr.size();

        if(limit > 50) {
            limit = 50;
        }

        List<Double> results = new ArrayList<>();

        if(genotype_curr.size()==0) {
            results.add(diversity);
            results.add(tmrca1);
            //results.add(geneticDiversity);
            if(getDivergence) {
                results.add(divergence);
            }

            return results;
        }


        for (int i = 0; i < limit; i++) {

            Integer vA = genotype_curr.get(Uniform.staticNextIntFromTo(0, genotype_curr.size()-1));
            Integer vB = genotype_curr.get(Uniform.staticNextIntFromTo(0, genotype_curr.size()-1));

            //if (vA.intValue() != (int)Double.NEGATIVE_INFINITY && vB.intValue() != (int)Double.NEGATIVE_INFINITY) {
            double dist = distance(vA, vB);
            //double geneticDist = geneticDistance(vA, vB);
            double div_a = patch_history.getBirth(vA);
            double div_b = patch_history.getBirth(vB);

            if(Double.isInfinite(div_a)) {
                div_a = 0.0;
            }
            if(Double.isInfinite(div_b)) {
                div_b = 0.0;
            }

            diversity += dist;
            //geneticDiversity += geneticDist;
            divergence += div_a;
            divergence += div_b;


            if (dist > tmrca1) {
                tmrca1 = dist;
            }
            sampleCount1 += 1;
        }

        if (sampleCount1 > 0) {
            diversity /= (double) sampleCount1;
            //geneticDiversity /= (double) sampleCount1;
            divergence /= (double) sampleCount1*2;



        }

        tmrca1 /= 2.0;

        results.add(diversity);
        results.add(tmrca1);
        //results.add(geneticDiversity);
        if(getDivergence) {
            results.add(divergence);
        }


        return results;

    }

    public byte[] serializeData(DefaultCoder coder) {


        ArrayList<Object> list = new ArrayList<>();

        list.add(genotype_curr_per_patch);
        list.add(curr_in_body);
        list.add(occupyPatchList);
        list.add(emptyPatchList);
        list.add(X_curr_i);
        list.add(Y_curr_i);
        list.add(patch_history.genotype);
        list.add(patch_history.mutationsFromParent);
        list.add(patch_history.birth);
        list.add(patch_history.death);
        list.add(patch_history.parent);
        list.add(patch_history.prevalence);
        list.add(patch_history.patch);
        list.add(patch_history.beta);
        list.add(genotype);
        list.add(timeThreshold_curr_i);
        list.add(timeRefactory_curr_i);



        byte serialized[] = coder.toByteArray(list);
//            Object deserialized = coder.toObject( serialized );
//
//            ByteOutputStream out = new ByteOutputStream();
//            bytes = conf.asByteArray(list);
//            FSTObjectOutput stream = new FSTObjectOutput(out);
//            stream.writeInt(bytes.length);
//            stream.write(bytes);
//            stream.flush();
//            out.flush();
//            stream.close();
//            out.close();
//            //FileOutputStream fileOut = new FileOutputStream("history.ser");
//            ByteArrayOutputStream bos = new ByteArrayOutputStream();
//            //bos.writeTo(fileOut);
//            ObjectOutputStream out = new ObjectOutputStream(bos);
//            out.writeObject(list);
//            out.flush();
//            bytes = bos.toByteArray();
//            out.close();
//            bos.close();
        //System.out.println("Serialized data is saved in bytes "+bytes.length);
        return serialized;


    }

    public void deserializaData(byte[] serialized, DefaultCoder coder) {

//        System.out.println(">"+bytes[0]);
        //ArrayList<Object> deserialized = null;
        // Deserialize in to new class object
//        Object deserialized = null;
//        ByteArrayInputStream bis = new ByteArrayInputStream(bytes);
//        ObjectInput in = null;


        Object readObj = coder.toObject(serialized);

//        try {
//
//            while (true) {
//                try {
//
//                    ByteArrayInputStream stream = new ByteArrayInputStream(bytes);
//                    FSTObjectInput in = new FSTObjectInput(stream);
////                    int len = in.readInt();
////                    int orglen = len;
////                    byte[] buffer = new byte[len]; // this could be reused !
////                    while (len > 0) {
////                        len -= in.read(buffer, buffer.length - len, len);
////                    }
//                    readObj = in.readObject();
//
//                    in.close();
// skipped: check for stream close

        ArrayList<Object> x = (ArrayList<Object>) readObj;
        this.genotype_curr_per_patch = (HashMap<Integer, List<Integer>>) x.get(0);
        this.curr_in_body = (ArrayList<Integer>) x.get(1);
        this.occupyPatchList = (ArrayList<Integer>) x.get(2);
        this.emptyPatchList = (ArrayList<Integer>) x.get(3);
        this.X_curr_i = (ArrayList<Integer>) x.get(4);
        this.Y_curr_i = (ArrayList<Integer>) x.get(5);
        this.patch_history.genotype = ((ArrayList<Integer>) x.get(6));
        this.patch_history.mutationsFromParent = ((ArrayList<Integer>) x.get(7));

        this.patch_history.birth = ((ArrayList<Double>) x.get(8));
        this.patch_history.death = ((ArrayList<Double>) x.get(9));
        this.patch_history.parent = ((ArrayList<Integer>) x.get(10));
        this.patch_history.prevalence = ((ArrayList<List<Integer>>) x.get(11));
        this.patch_history.patch = ((ArrayList<Integer>) x.get(12));
        this.patch_history.beta = ((ArrayList<Double>) x.get(13));
        this.genotype = (Integer) x.get(14);
        this.timeThreshold_curr_i = ((ArrayList<Double>) x.get(15));
        this.timeRefactory_curr_i = ((ArrayList<Double>) x.get(16));
//                    bis.close();
//                    in.close();
//                } catch (EOFException ex1) {
//                    break;
//                }
//            }
//        } catch (IOException e) {
//            e.printStackTrace();
//        }catch (ClassNotFoundException e) {
//            e.printStackTrace();
//        }
    }

    public static void main(String [] args) throws  MPIException{

        patchSimMPI test = new patchSimMPI();
        params inputParams = new params();
        inputParams.runTime = Double.parseDouble(args[0]);
        inputParams.startTime = Double.parseDouble(args[1]);
        inputParams.Npatches = Integer.parseInt(args[2]);
        inputParams.S = Integer.parseInt(args[3]);
        inputParams.I = Integer.parseInt(args[4]);
        inputParams.nu = Double.parseDouble(args[5]);
        inputParams.c = Double.parseDouble(args[6]);
//        inputParams.U = Double.parseDouble(args[7]);
        MPI.Init(args);
        test.runSim(inputParams);
        MPI.Finalize();


    }

    private int sumArray(int[] array) {

        int i = 0;
        int sum = 0;
        while(i < array.length) {

            sum += array[i];
            i++;
        }
        return sum;
    }

    private  Comparator<Integer>  comparator1 = new Comparator<Integer>() {
        @Override
        public int compare(Integer o1, Integer o2) {


            List<Integer> prevalence1 = patch_history.prevalence.get(o1);
            List<Integer> prevalence2 = patch_history.prevalence.get(o2);

            double diff = patch_history.getBeta(o1) - patch_history.getBeta(o2);
            if(diff == 0) {

                return 0;
            }
            else if(diff > 0) {
                return 1;
            }
            else{
                return -1;
            }
        }
    };

    public int chooseBin(double [] prob, double randomNo) {

        int bin = 0;

        while(prob[bin] < randomNo) {

            bin++;

            if(bin==prob.length){
                break;
            }

        }

        return bin;

    }

    private double getRandomNo() {

        Uniform randomNo = new Uniform(0,1, params.randomGenerator);

        return randomNo.nextDouble();

    }
}
