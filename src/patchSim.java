/**
 * Created by jayna on 17/05/2018.
 */

import cern.colt.list.DoubleArrayList;
import cern.jet.random.Exponential;
import cern.jet.random.Poisson;
import cern.jet.random.Uniform;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import mpi.*;

import javax.xml.crypto.Data;


public class patchSim {

    params params;
    infectionHistory patch_history;

    double tau; // time step
    int genotype;
    Map<Integer, List<Integer>> genotype_curr_per_patch;
    List<Double> potentialSamplingTimepoints;


    public patchSim() {

        tau = 0.25;
        patch_history = new infectionHistory();

        genotype = 1;
        genotype_curr_per_patch = new HashMap<>();
        //potentialSamplingTimepoints = new ArrayList<>();

    }


    public void runSim(params params) {

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
            writer1.write("Time\tPatch_no\tTotal_infected\tUnique_genotypes\tGenealogical_diversity\tTMRCA\tGenetic_diversity\tGenetic_divergence\tSim\n");

            writer2 = new FileWriter(new File(summaryFile));
            writer2.write("Time\tTotal_infected\tUnique_genotypes\tOccupied_patches\tGenealogical_diversity\tTMRCA\tGenetic_diversity\tGenetic_divergence\tSim\n");
        } catch (IOException e) {
            e.printStackTrace();
        }

        double t_curr = tau;
        double t_max = params.runTime;
        int Y_curr = params.I;
        int X_curr = params.S;
        DoubleArrayList withinPatchRates = new DoubleArrayList();
        DoubleArrayList betweenPatchRates = new DoubleArrayList();
        List<Integer> curr_in_body = new ArrayList<>();


        int[] Y_curr_i = new int[params.Npatches];
        int[] X_curr_i = new int[params.Npatches];
        double[] beta_i = new double[params.Npatches];
        double[] timeThreshold_curr_i = new double[params.Npatches];
        double[] timeRefactory_curr_i = new double[params.Npatches];
        int[] daysRefactory_i = new int[params.Npatches];
        int[] threshold_i = new int[params.Npatches];
        double[] colonize_i = new double[params.Npatches];
        double[] extinct_i = new double[params.Npatches];

        List<Integer> emptyPatchList = new ArrayList<>();
        List<Integer> occupyPatchList = new ArrayList<>();



        //initializing patches with number of infected and susceptible cells

        Arrays.fill(Y_curr_i, 0);
        Arrays.fill(X_curr_i, (X_curr+Y_curr));
        Arrays.fill(timeThreshold_curr_i, 0.0);
        Arrays.fill(timeRefactory_curr_i, 0);
        Arrays.fill(beta_i, params.beta);
        Arrays.fill(daysRefactory_i, 5);





        double colonize = 0.0;
        double extinct = 0.0;
        double death = 0.1;


        //int i1 = (int)Math.ceil(Uniform.staticNextDouble()*(double)Y_curr_i.length);
//        int i2 = (int)Math.ceil(Uniform.staticNextDouble()*(double)Y_curr_i.length);

        Y_curr_i[0] += Y_curr;
        X_curr_i[0] -= Y_curr;


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
            //daysRefactory_i[i] = Poisson.staticNextInt(5);
            //beta_i[i] = Exponential.staticNextDouble(params.beta);//Uniform.staticNextDoubleFromTo((0.1*params.beta),(1.9*params.beta));


        }


        for (int t = 1; t < (t_max / tau); t++) {

            t_curr += tau;

            int total_infected = 0;
            int total_genotypes = 0;

            colonize = params.c;
            extinct = params.nu;
            death = 0.1;

            if(t_curr < 40) {

                extinct = 0;
                colonize = 0.1;

            }

            String sim = "ext_"+ extinct+"_col_"+colonize+"_Npatch_"+params.Npatches+"_PatchSize_"+(params.S+params.I);

            for(Integer patch: genotype_curr_per_patch.keySet()) {

                total_infected+=genotype_curr_per_patch.get(patch).size();
                total_genotypes+=(new HashSet<>(genotype_curr_per_patch.get(patch))).size();
                if(genotype_curr_per_patch.get(patch).size()==0) {
                    emptyPatchList.add(patch);
                }
                else{
                    occupyPatchList.add(patch);
                }



                if(t_curr % 10 == 0) {
                    List<Double> diversity_results = updateDiversity(genotype_curr_per_patch.get(patch), false);

                    double diversity = diversity_results.get(0);
                    double tmrca = diversity_results.get(1);
                    double geneticDiversity = diversity_results.get(2);

                    try {
                        writer1.write(t_curr + "\t" + (patch + 1) + "\t" + genotype_curr_per_patch.get(patch).size() +
                                "\t" + (new HashSet<>(genotype_curr_per_patch.get(patch))).size() +
                                "\t" + diversity + "\t" + tmrca + "\t" + geneticDiversity + "\t" + sim + "\n");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }

            }


            for (int i = 0; i < patch_history.prevalence.size(); i++) {

                //List<Integer> prevalence = patch_history.prevalence.get(i);
                int prev_curr = patch_history.prevalence.get(i).get(t - 1);
                patch_history.prevalence.get(i).set(t, prev_curr);


                Integer genotype_i = patch_history.genotype.get(i);

                if (patch_history.prevalence.get(i).get(t) > 0) {

                    Integer[] genotype_array = new Integer[patch_history.prevalence.get(i).get(t - 1)];
                    Arrays.fill(genotype_array, genotype_i);
                    curr_in_body.addAll(Arrays.asList(genotype_array));
                }

            }


            if(t_curr%10==0) {
                List<Double> globalDiversity = updateDiversity(curr_in_body, true);
                double diversity = globalDiversity.get(0);
                double tmrca = globalDiversity.get(1);
                double geneticDiversity = globalDiversity.get(2);
                double divergence = globalDiversity.get(3);

                System.out.println("t_curr \t" + t_curr + "\ttotal_infected: \t" + Math.log10(total_infected) + "\ttotal_genotypes: \t" + total_genotypes +
                        "\tpatch_occupied: \t"+ occupyPatchList.size()+ "\t"+geneticDiversity + "\t"+divergence);



                writeOutput(writer2, t_curr, total_infected, total_genotypes, occupyPatchList.size(), globalDiversity, sim);

            }

            betweenPatchRates.add(extinct*occupyPatchList.size());
            betweenPatchRates.add(colonize*emptyPatchList.size());


            Poisson poisson;
            int j;
            for(int event = 0; event < betweenPatchRates.size(); event++) {

                poisson = new Poisson(tau * betweenPatchRates.get(event), params.randomGenerator);
                int num = poisson.nextInt();

                if(event == 0) {

                    j = 0;

                    //extinction

//                    List<Integer> occupyPatchList_copy = deepCopy(occupyPatchList);
//                    List<Integer> emptyPatchList_copy = deepCopy(emptyPatchList);


                    int min = Math.min(num, occupyPatchList.size());

                    while (j < min) {

                        // patch extinction

                        int r = Uniform.staticNextIntFromTo(0, occupyPatchList.size()-1);
                        int patch_index = occupyPatchList.get(r);

                        double prob_extinct = 0;

                        if(Y_curr_i[patch_index] >= threshold_i[patch_index]) {

//                            prob_extinct = Uniform.staticNextDouble();
//                        }
//
//                        if(prob_extinct > 0.5) {

                            Y_curr_i[patch_index] = 0;
                            X_curr_i[patch_index] = params.S + params.I;

                            //randomly choose a genotype to seed/colonize new patch

                            timeThreshold_curr_i[patch_index] = 0.0;
                            timeRefactory_curr_i[patch_index] = t_curr;

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


                            occupyPatchList.remove((Integer)patch_index);
                            emptyPatchList.add(patch_index);

                        }

                        j++;
                    }

//                    occupyPatchList = deepCopy(occupyPatchList_copy);
//                    emptyPatchList = deepCopy(emptyPatchList_copy);


                }
                else{

                    // patch colonisation

                    j = 0;

//                    List<Integer> emptyPatchList_copy = deepCopy(emptyPatchList);
//                    List<Integer> occupyPatchList_copy = deepCopy(occupyPatchList);


                    int min = Math.min(num, emptyPatchList.size());

                    while (j < min) {

                        int r = Uniform.staticNextIntFromTo(0, emptyPatchList.size()-1);

                        int patch_index = emptyPatchList.get(r);

                        int probTimeRefactory = daysRefactory_i[patch_index]; // choose probabilistically how many days patch remains in refactory period (make a variable of the model)

                        if((t_curr-timeRefactory_curr_i[patch_index]) > (double)probTimeRefactory) {

                            if(curr_in_body.size() == 0) {

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


                            double rand = Uniform.staticNextDouble();

                            double prob_colonize = 0.5;

//
//                            if (rand < prob_colonize) {

                            timeRefactory_curr_i[patch_index] = 0.0;

                            int newInfections = Uniform.staticNextIntFromTo(1, 3);

                            Y_curr_i[patch_index] = newInfections;
                            X_curr_i[patch_index] -= newInfections;

                            r = Uniform.staticNextIntFromTo(0, curr_in_body.size()-1);
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
                            //curr_in_body.clear();

                            emptyPatchList.remove((Integer)patch_index);
                            occupyPatchList.add((Integer)(patch_index));

                        }
                        j++;
                    }


//                    occupyPatchList = deepCopy(occupyPatchList_copy);
//                    emptyPatchList = deepCopy(emptyPatchList_copy);

                }
            }
            betweenPatchRates.removeAll(betweenPatchRates);


            for (int p = 0; p < params.Npatches; p++) {

                int X_temp = X_curr_i[p];
                int Y_temp = Y_curr_i[p];


                //withinPatchRates.add(mutate * Y_temp/total); // mutation rate
                withinPatchRates.add((beta_i[p] * X_temp * Y_temp) / (X_temp+Y_temp)); // infection rate
//                withinPatchRates.add((params.nu*(Math.pow(Y_temp,5)/Math.pow(threshold_i[p], 5)))); // extinction
//                withinPatchRates.add((((double)total-(double)Y_temp)/total)*((2.0-(timeRefactory_curr_i[p]/t_curr)) * params.c)); // colonise
                withinPatchRates.add(death*(Y_temp));


                //Poisson poisson;
                int noOfRates = withinPatchRates.size();
                //int j;


                //events in patches that affect X and Y individuals
                for (int event = 0; event < noOfRates; event++) {

                    poisson = new Poisson(tau * withinPatchRates.get(event), params.randomGenerator);
                    int num = poisson.nextInt();

                    int patch_index = p;

                    List<Integer> genotype_curr = new ArrayList<>();
                    genotype_curr.addAll(genotype_curr_per_patch.get(patch_index));


                    if (event == 0) {
                        // mutations in patch_i

//                        j = 0;
//
//                        int min = Math.min(num, Y_curr_i[patch_index]);
//
//                        while (j < min) {
//
//                            //choose genotype
//                            int index = (int) Math.floor(Uniform.staticNextDouble() * genotype_curr.size());
//                            Integer chosen_genotype = genotype_curr.get(index);
//                            int parent_index = patch_history.genotype.indexOf(chosen_genotype);
//
//
//                            patch_history.logParent(chosen_genotype); // do we want to log parent index or genotype?
//                            patch_history.logGenotype(genotype);
//                            patch_history.logBirth(t_curr);
//                            patch_history.logDeath(Double.POSITIVE_INFINITY);
//                            patch_history.logPatch(patch_index);
//                            patch_history.logParentOrigin(0);
//
//                            List<Integer> prevalence = new ArrayList<>();
//                            prevalence.addAll(patch_history.prevalence.get(parent_index));
//
//                            //prevalence is bit off need to check this
//                            int prevalence_curr = prevalence.get(t - 1);
//                            int prevalence_next = prevalence.get(t);
//
//                            if (prevalence_next > 0) {
//
//                                prevalence_curr = prevalence_next;
//                            }
//
//                            if((prevalence_curr-1)==0) {
//
//                                patch_history.death.set(parent_index, t_curr);
//                            }
//                            patch_history.prevalence.get(parent_index).set(t, prevalence_curr - 1);
//
//                            List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[prevalence.size()]));
//                            Collections.fill(new_prevalence, 0);
//                            new_prevalence.set(t, 1);
//                            patch_history.prevalence.add(new_prevalence);
//
//                            genotype_curr.remove(chosen_genotype);
//                            genotype_curr.add(genotype);
//                            genotype++;
//
//                            j++;
//                        }
//
//
//                    }
//                    else { // new infection


                        int min = Math.min(X_curr_i[patch_index], num);

                        j = 0;


                        while (j < min) {

                            int mutations = Poisson.staticNextInt(params.U); //double check if rate vs mean


                            if(mutations == 0) {
                                int r = (int) Math.floor(Uniform.staticNextDouble() * genotype_curr.size());
                                Integer chosen_genotype = genotype_curr.get(r);
                                int genotype_index = patch_history.genotype.indexOf(chosen_genotype);

//                                List<Integer> prevalence = new ArrayList<>();
//                                prevalence.addAll(patch_history.prevalence.get(genotype_index));

                                int prevalence_curr = patch_history.prevalence.get(genotype_index).get(t - 1);
                                int prevalence_next = patch_history.prevalence.get(genotype_index).get(t);

                                X_curr_i[patch_index] = X_curr_i[patch_index] - 1;
                                Y_curr_i[patch_index] = Y_curr_i[patch_index] + 1;

                                if (Y_curr_i[patch_index] >= threshold_i[patch_index]) {

                                    timeThreshold_curr_i[patch_index] = t_curr;

                                }

                                if (prevalence_next > 0) {

                                    prevalence_curr = prevalence_next;
                                }

                                genotype_curr.add(chosen_genotype);
                                patch_history.prevalence.get(genotype_index).set(t, prevalence_curr + 1);
                            }
                            else {


                                //choose genotype
                                int index = (int) Math.floor(Uniform.staticNextDouble() * genotype_curr.size());
                                Integer chosen_genotype = genotype_curr.get(index);
                                int parent_index = patch_history.genotype.indexOf(chosen_genotype);


                                patch_history.logParent(chosen_genotype); // do we want to log parent index or genotype?
                                patch_history.logGenotype(genotype);
                                patch_history.logMutations(mutations);
                                patch_history.logBirth(t_curr);
                                patch_history.logDeath(Double.POSITIVE_INFINITY);
                                patch_history.logPatch(patch_index);

//                                List<Integer> prevalence = new ArrayList<>();
//                                prevalence.addAll(patch_history.prevalence.get(parent_index));

                                //prevalence is bit off need to check this
                                int prevalence_curr = patch_history.prevalence.get(parent_index).get(t - 1);
                                int prevalence_next = patch_history.prevalence.get(parent_index).get(t);

                                if (prevalence_next > 0) {

                                    prevalence_curr = prevalence_next;
                                }

                                if((prevalence_curr-1)==0) {

                                    patch_history.death.set(parent_index, t_curr);
                                }
                                patch_history.prevalence.get(parent_index).set(t, prevalence_curr - 1);

                                List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[patch_history.prevalence.get(parent_index).size()]));
                                Collections.fill(new_prevalence, 0);
                                new_prevalence.set(t, 1);
                                patch_history.prevalence.add(new_prevalence);

                                genotype_curr.remove(chosen_genotype);
                                genotype_curr.add(genotype);
                                genotype++;

                            }
                            j++;
                        }

                    }
                    else{

                        int min = Math.min(Y_curr_i[patch_index], num);

                        j = 0;

                        while(j < min) {

                            X_curr_i[patch_index] = X_curr_i[patch_index] + 1;
                            Y_curr_i[patch_index] = Y_curr_i[patch_index] - 1;

                            int index = (int) Math.floor(Uniform.staticNextDouble() * genotype_curr.size());
                            Integer chosen_genotype = genotype_curr.get(index);
                            int genotype_index = patch_history.genotype.indexOf(chosen_genotype);

//                            List<Integer> prevalence = new ArrayList<>();
//                            prevalence.addAll(patch_history.prevalence.get(genotype_index));

                            int prevalence_curr = patch_history.prevalence.get(genotype_index).get(t);
//                            int prevalence_next = prevalence.get(t);
//
//                            if (prevalence_next > 0) {
//
//                                prevalence_curr = prevalence_next;
//                            }

                            if((prevalence_curr-1)==0) {

                                patch_history.death.set(genotype_index, t_curr);
                            }

                            patch_history.prevalence.get(genotype_index).set(t, prevalence_curr - 1);
                            genotype_curr.remove(chosen_genotype);
                            j++;


                        }



                    }


                    genotype_curr_per_patch.replace(patch_index, genotype_curr);

                }

                withinPatchRates.removeAll(withinPatchRates);

//
//                int size = 0;
//                for (int pr = 0; pr < patch_history.prevalence.size(); pr++) {
//
//                    if (patch_history.patch.get(pr) == p) {
//                        size += patch_history.prevalence.get(pr).get(t);
//                    }
//                }

            }

            curr_in_body.clear();
            occupyPatchList.clear();
            emptyPatchList.clear();

        }

        try {
            writer1.close();
            writer2.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }



    private void initialiseHistory(infectionHistory history, int [] Y_curr_i, double t_max, int n_patches) {

        int i = 0;
        while(i < n_patches) {

            if(Y_curr_i[i] > 0) {
                history.logGenotype(genotype);
                history.logMutations(0);
                history.logBirth(Double.NEGATIVE_INFINITY);
                history.logDeath(Double.POSITIVE_INFINITY);
                history.logParent((int) Double.NEGATIVE_INFINITY);
                List<Integer> prevalence_through_time = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau)]));
                Collections.fill(prevalence_through_time, 0);
                prevalence_through_time.set(0, Y_curr_i[i]);
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
            writer.write(t+"\t"+total_infected+"\t"+genotypes+"\t"+n_patches+"\t"+diversity.get(0)+"\t"+diversity.get(1) + "\t"+diversity.get(2)+ "\t"+diversity.get(3) + "\t"+sim+"\n");
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
//        }
//        else{
//            return 0;
//
//        }
    }



    public List<Double> updateDiversity(List<Integer> genotype_curr, boolean getDivergence) {

        double diversity = 0.0;
        double tmrca1 = 0.0;
        int sampleCount1 = 0;

        double geneticDiversity = 0.0;
        double divergence = 0.0;

//        List<Integer> indices1 = new ArrayList<>();
//        List<Integer> indices2 = new ArrayList<>();


        int limit = genotype_curr.size();

        if(limit > 50) {
            limit = 50;
        }

//        for(int i= 0; i < genotype_curr.size(); i++) {
//
//            indices1.add(i);
//            indices2.add(i);
//        }

//        Collections.shuffle(indices1, new Random((long)(Math.floor(Math.random()*10000))));
//        Collections.shuffle(indices2, new Random((long)(Math.floor(Math.random()*10000))));
        List<Double> results = new ArrayList<>();

        if(genotype_curr.size()==0) {
            results.add(diversity);
            results.add(tmrca1);
            results.add(geneticDiversity);
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
            double geneticDist = geneticDistance(vA, vB);
            double div_a = patch_history.getMutationsFromParent(vA, 1);
            double div_b = patch_history.getMutationsFromParent(vB, 1);

            diversity += dist;
            geneticDiversity += geneticDist;
            divergence += div_a;
            divergence +=div_b;


            if (dist > tmrca1) {
                tmrca1 = dist;
            }
            sampleCount1 += 1;
        }
        //}
        if (sampleCount1 > 0) {
            diversity /= (double) sampleCount1;
            geneticDiversity /= (double) sampleCount1;
            divergence /= (double) sampleCount1*2;



        }

        tmrca1 /= 2.0;

        results.add(diversity);
        results.add(tmrca1);
        results.add(geneticDiversity);
        if(getDivergence) {
            results.add(divergence);
        }


        return results;

    }

    public static void main(String [] args) throws MPIException {

//        patchSim test = new patchSim();
//        params inputParams = new params();
//        inputParams.runTime = Double.parseDouble(args[0]);
//        inputParams.startTime = Double.parseDouble(args[1]);
//        inputParams.Npatches = Integer.parseInt(args[2]);
//        inputParams.S = Integer.parseInt(args[3]);
//        inputParams.I = Integer.parseInt(args[4]);
//        inputParams.nu = Double.parseDouble(args[5]);
//        inputParams.c = Double.parseDouble(args[6]);
//        inputParams.U = Double.parseDouble(args[7]);
//        test.runSim(inputParams);


// Hello World
//        // number of processes
//        int size;
//        // process rank
//        int rank;
//
//        // initialise MPI
//        MPI.Init(args) ;
//        // retrieve size of default communicator
//        size = MPI.COMM_WORLD.getSize();
//        // retrieve process rank
//        rank = MPI.COMM_WORLD.getRank();
//        // print message
//        System.out.println("Hello world from process "+ rank+ " of "+ size);
//        // finalise MPI
//        MPI.Finalize();



        MPI.Init(args);

        int rank = MPI.COMM_WORLD.getRank() ; //The current process.
        int size = MPI.COMM_WORLD.getSize() ; //Total number of processes
        int peer ;

        params[] buffer  = new params[10];
        params x = new params();
        int len = 1 ;
        int dataToBeSent = 99 ;
        int tag = 100 ;

        if(rank == 0) {

            buffer[0] = x ;
            peer = 1 ;
            MPI.COMM_WORLD.send(buffer, len, MPI.INT, peer, tag) ;
            System.out.println("process <"+rank+"> sent a msg to process <"+peer+">") ;

        } else if(rank == 1) {

            peer = 0 ;
            Status status = MPI.COMM_WORLD.recv(buffer, buffer.length, MPI.INT, peer, tag);
            System.out.println("process <"+rank+"> recv'ed a msg\n" +"\tdata   <"+buffer[0].Npatches    +"> \n"+
                    "\tsource <"+status.getSource()+"> \n"+"\ttag    <"+status.getTag()   +"> \n"+"\tcount  <"+status.getCount(MPI.INT) +">") ;

        }

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

            double diff = patch_history.getFitness(o1) - patch_history.getFitness(o2);
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
