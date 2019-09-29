/**
 * Created by jayna on 17/05/2018.
 */

import cern.colt.list.DoubleArrayList;
import cern.jet.random.Exponential;
import cern.jet.random.Gamma;
import cern.jet.random.Poisson;
import cern.jet.random.Uniform;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.nio.Buffer;
import java.util.*;
import mpi.*;
import org.apache.commons.lang3.ArrayUtils;

import java.io.*;
import java.nio.ByteBuffer;
import java.util.stream.DoubleStream;


public class patchSim {

    params params;
    infectionHistory patch_history;

    double tau; // time step
    int genotype;
    Map<Integer, List<Integer>> genotype_curr_per_patch;
    List<Double> potentialSamplingTimepoints;
    List<Integer> emptyPatchList;
    List<Integer> occupyPatchList;
    List<Integer> thresholdList;

    int total_infected;
    int total_genotypes;
    int [] threshold_i;

    public patchSim() {

        tau = 0.5;
        patch_history = new infectionHistory();
        genotype = 1;
        genotype_curr_per_patch = new HashMap<>();
        emptyPatchList = new ArrayList<>();
        occupyPatchList = new ArrayList<>();
        thresholdList = new ArrayList<>();

        total_genotypes = 0;
        total_infected = 0;
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
        List<Integer> curr_in_body = new ArrayList<>();


        int[] Y_curr_i = new int[params.Npatches];
        //int[] X_curr_i = new int[params.Npatches];
        double[] beta_i = new double[params.Npatches];
        double[] timeThreshold_curr_i = new double[params.Npatches];
        double[] timeRefactory_curr_i = new double[params.Npatches];
        int[] daysRefactory_i = new int[params.Npatches];
        threshold_i = new int[params.Npatches];
        double[] colonize_i = new double[params.Npatches];
        double[] extinct_i = new double[params.Npatches];

//        List<Integer> emptyPatchList = new ArrayList<>();
//        List<Integer> occupyPatchList = new ArrayList<>();

        //initializing patches with number of infected and susceptible cells

        Arrays.fill(Y_curr_i, 0);
        //Arrays.fill(X_curr_i, (X_curr+Y_curr));
        Arrays.fill(timeThreshold_curr_i, 0.0);
        Arrays.fill(timeRefactory_curr_i, 0);
        //Arrays.fill(beta_i, params.beta);
        Arrays.fill(daysRefactory_i, 5);

        double colonize = 0.0;
        double extinct = 0.0;
        double death = 0.1;


        Y_curr_i[0] += Y_curr;
        //X_curr_i[0] -= Y_curr;

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

        int N = params.S+params.I;
        Arrays.fill(threshold_i, (int)Math.ceil(0.8*N));
        for(int i = 0; i < threshold_i.length; i++) {


            beta_i[i] = params.beta*Math.abs(Gamma.staticNextDouble(0.2, 1)-1);
            //threshold_i[i] = Uniform.staticNextIntFromTo((int)Math.ceil(0.3*N),(int)Math.ceil(N));
            if(genotype_curr_per_patch.get(i).size() == 0) {
                emptyPatchList.add(i);
            }
            else{
                occupyPatchList.add(i);
            }
        }

        System.out.println(DoubleStream.of(beta_i).min().toString());
        System.out.println(DoubleStream.of(beta_i).max().toString());




        for (int t = 1; t < (t_max / tau); t++) {

            t_curr += tau;


            colonize = params.c;
            extinct = params.nu;
            death = 0.1;

            if(t_curr < 40) {

                extinct = 0;
                //colonize = 0.1;

            }


            String sim = "ext_"+ extinct+"_col_"+colonize+"_Npatch_"+params.Npatches+"_PatchSize_"+(params.S+params.I);


            updatePatchStatistics(writer1, t_curr, sim);

            //determine what is present currently across all patches
            update_curr_in_body(curr_in_body, t);


            if(t_curr%10==0) {
                List<Double> globalDiversity = updateDiversity(curr_in_body, true);
                double diversity = globalDiversity.get(0);
                double tmrca = globalDiversity.get(1);
                //double geneticDiversity = globalDiversity.get(2);
                double divergence = globalDiversity.get(2);

                System.out.println("t_curr \t" + t_curr + "\ttotal_infected: \t" + Math.log10(total_infected) + "\ttotal_genotypes: \t" + total_genotypes +
                        "\tpatch_occupied: \t"+ occupyPatchList.size() + "\t"+(tmrca));




                writeOutput(writer2, t_curr, total_infected, total_genotypes, occupyPatchList.size(), globalDiversity, sim);

            }

            double N_occupied = (double)occupyPatchList.size();


            betweenPatchRates.add(extinct*thresholdList.size());
            betweenPatchRates.add(colonize*(params.Npatches-N_occupied));


            Poisson poisson;
            int j;
            for(int event = 0; event < betweenPatchRates.size(); event++) {

                poisson = new Poisson(tau * betweenPatchRates.get(event), params.randomGenerator);
                int num = poisson.nextInt();

                if(event == 0) {

                    j = 0;

                    //extinction

                    int min = Math.min(num, thresholdList.size());

                    while (j < min) {

                        // patch extinction

                        int r = Uniform.staticNextIntFromTo(0, thresholdList.size()-1);
                        int patch_index = thresholdList.get(r);

                        Y_curr_i[patch_index] = 0;

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

                        genotype_curr_per_patch.get(patch_index).clear();
                        //occupyPatchList.remove((Integer)patch_index);
                        //thresholdList.remove((Integer)patch_index);
                        emptyPatchList.add(patch_index);


                        //}
                        j++;

                    }
                    thresholdList.removeAll(emptyPatchList);
                    occupyPatchList.removeAll(emptyPatchList);


                }
                else{



                    // patch colonisation

                    j = 0;

                    int min = Math.min(num, emptyPatchList.size());

                    if(min == 0) {
                        break;
                    }

                    while (j < min) {

                        int r = Uniform.staticNextIntFromTo(0, emptyPatchList.size()-1);

                        int patch_index = emptyPatchList.get(r);

                        //int probTimeRefactory = daysRefactory_i[patch_index]; // choose probabilistically how many days patch remains in refactory period (make a variable of the model)

                        //if((t_curr-timeRefactory_curr_i[patch_index]) > (double)probTimeRefactory) {

                        if(curr_in_body.size() == 0) {

                            break;
                        }

//                        int s = Uniform.staticNextIntFromTo(0, occupyPatchList.size()-1);
//                        int source = occupyPatchList.get(s);

                        // Negative frequency dependent selection of colonizing patch

                        Collections.sort(occupyPatchList, comparator1);
                        double[] prob = new double[occupyPatchList.size()];
                        double[] relProb = new double[occupyPatchList.size()];

                        int a = 0;
                        for(Integer o: occupyPatchList) {

                            double fitness = (double)(total_infected-Y_curr_i[o])/(double)total_infected;
                            prob[a] = fitness;

                            if(a > 0) {
                                prob[a] = prob[a] + prob[a-1];
                            }
                            a++;

                        }

                        a = 0;

                        double max = prob[prob.length-1];
                        for(double d: prob) {

                            relProb[a] = d/max;
                            a++;

                        }

                        int index = 0;
                        double freq = relProb[index];
                        double rand = Uniform.staticNextDouble();
                        while(freq < rand) {
                            if (index == (relProb.length - 1)) {
                                break;
                            }
                            index++;
                            freq = relProb[index];
                        }

                        int source = occupyPatchList.get(index);

                        r = Uniform.staticNextIntFromTo(0, genotype_curr_per_patch.get(source).size()-1);

                        //Collections.shuffle(curr_in_body);

                        //timeRefactory_curr_i[patch_index] = 0.0;

                        int newInfections = Uniform.staticNextIntFromTo(1, 3);

                        Y_curr_i[patch_index] = newInfections;
                        //X_curr_i[patch_index] -= newInfections;


                        //r = Uniform.staticNextIntFromTo(0, curr_in_body.size()-1);

                        Integer colonizing_genotype = genotype_curr_per_patch.get(source).get(r);  //curr_in_body.get(r);

                        for (int i = 0; i < newInfections; i++) {

                            genotype_curr_per_patch.get(patch_index).add(genotype);

                        }
                        patch_history.logBirth(t_curr);
                        patch_history.logDeath(Double.POSITIVE_INFINITY);
                        patch_history.logGenotype(genotype);
                        //patch_history.logFitness(1.0-(double)Y_curr_i[patch_index]/(double)N);
                        //patch_history.logMutations(0);
                        patch_history.logParent(colonizing_genotype);
                        patch_history.logPatch(patch_index);

                        List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau)]));
                        Collections.fill(new_prevalence, 0);
                        new_prevalence.set(t, newInfections);
                        patch_history.prevalence.add(new_prevalence);
                        genotype++;

                        emptyPatchList.remove((Integer)patch_index);
                        occupyPatchList.add(patch_index);

                        //tempOccupyList.add(patch_index);
                        //}
                        j++;
                    }
                    //emptyPatchList.removeAll(tempOccupyList);
                    //occupyPatchList.addAll(tempOccupyList);



                }
            }
            betweenPatchRates.removeAll(betweenPatchRates);


            for (int p = 0; p < params.Npatches; p++) {

                //int X_temp = X_curr_i[p];
                int Y_temp = Y_curr_i[p];

                withinPatchRates.add((beta_i[p] * (N-Y_temp) * Y_temp) / (double)(N)); // infection rate
                withinPatchRates.add(death*(Y_temp));


                //Poisson poisson;
                int noOfRates = withinPatchRates.size();


                //events in patches that affect X and Y individuals
                for (int event = 0; event < noOfRates; event++) {

                    poisson = new Poisson(tau * withinPatchRates.get(event), params.randomGenerator);
                    int num = poisson.nextInt();

                    int patch_index = p;

                    List<Integer> genotype_curr = new ArrayList<>();
                    genotype_curr = (genotype_curr_per_patch.get(patch_index));


                    if (event == 0) {

                        // infection

                        int min = Math.min((N-Y_curr_i[patch_index]), num);

                        j = 0;


                        while (j < min) {

                            //determine number of mutations
                            int mutations = Poisson.staticNextInt(params.U); //double check if rate vs mean


                            if(mutations == 0) {
                                int r = (int) Math.floor(Uniform.staticNextDouble() * genotype_curr.size());
                                Integer chosen_genotype = genotype_curr.get(r);
                                int genotype_index = patch_history.genotype.indexOf(chosen_genotype);

//                                List<Integer> prevalence = new ArrayList<>();
//                                prevalence.addAll(patch_history.prevalence.get(genotype_index));

                                int prevalence_curr = patch_history.prevalence.get(genotype_index).get(t - 1);
                                int prevalence_next = patch_history.prevalence.get(genotype_index).get(t);

                                //X_curr_i[patch_index] = X_curr_i[patch_index] - 1;
                                Y_curr_i[patch_index] = Y_curr_i[patch_index] + 1;


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
                                //patch_history.logMutations(mutations);
                                //patch_history.logFitness(1.0-(double)Y_curr_i[patch_index]/(double)N);
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

                            if(Y_curr_i[patch_index] >= threshold_i[patch_index] && !thresholdList.contains(new Integer(patch_index))) {


                                thresholdList.add(patch_index);
                            }

                            j++;
                        }

                    }
                    else{

                        int min = Math.min(Y_curr_i[patch_index], num);

                        j = 0;

                        while(j < min) {

                            //X_curr_i[patch_index] = X_curr_i[patch_index] + 1;
                            Y_curr_i[patch_index] = Y_curr_i[patch_index] - 1;

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

                            if((prevalence_curr-1)==0) {

                                patch_history.death.set(genotype_index, t_curr);
                            }

                            patch_history.prevalence.get(genotype_index).set(t, prevalence_curr - 1);
                            genotype_curr.remove(chosen_genotype);

                            if(Y_curr_i[patch_index] == 0) {

                                occupyPatchList.remove((Integer)patch_index);
                                thresholdList.remove((Integer)patch_index);
                                emptyPatchList.add(patch_index);
                            }

                            j++;

                        }

                    }




                    genotype_curr_per_patch.replace(patch_index, genotype_curr);


                }

                withinPatchRates.removeAll(withinPatchRates);

            }

            curr_in_body.clear();
//            occupyPatchList.clear();
//            emptyPatchList.clear();

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
                //patch_history.logFitness(1-(double)Y_curr_i[i]/(double)(params.S+params.I));
                //history.logMutations(0);
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
//        }
//        else{
//            return 0;
//
//        }
    }

    public void updatePatchStatistics(FileWriter writer, double t_curr, String sim) {

        total_infected = 0;
        total_genotypes = 0;

        //for(Integer patch: genotype_curr_per_patch.keySet()) {


        for (Map.Entry<Integer, List<Integer>> entry : genotype_curr_per_patch.entrySet()) {

            Integer patch = entry.getKey();

            total_infected+=entry.getValue().size();
            total_genotypes+=(new HashSet<>(entry.getValue()).size());

            if(t_curr % 10 == 0) {
                List<Double> diversity_results = updateDiversity(genotype_curr_per_patch.get(patch), true);

                double diversity = diversity_results.get(0);
                double tmrca = diversity_results.get(1);
                double divergence = diversity_results.get(2);

                try {
                    writer.write(t_curr + "\t" + (patch + 1) + "\t" + genotype_curr_per_patch.get(patch).size() +
                            "\t" + (new HashSet<>(genotype_curr_per_patch.get(patch))).size() +
                            "\t" + diversity + "\t" + tmrca + "\t" + divergence + "\t" + sim + "\n");
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

                List<Integer> tmp = Collections.nCopies(prev_curr, genotype_i);
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

    public int choosePatch(double[] prob) {

        double x = Uniform.staticNextDouble();

        int index = 0;
        double i = prob[index];
        int patch = occupyPatchList.get(index);


        while(i < x || genotype_curr_per_patch.get(patch).size() == 0) {
            //System.out.println(parent_fitness+","+x);
            index++;
            if(index == (prob.length-1)) {
                break;
            }
            i = prob[index];
            patch = occupyPatchList.get(index);
            //x = Uniform.staticNextDouble();
        }

        //Double real_fitness = iMatrix.getFitness(parent_index);
        //System.out.println(patch+", "+genotype_curr_per_patch.get(patch).size());
        return occupyPatchList.get(index);
    }
    public static void main(String [] args) {

        patchSim test = new patchSim();
        params inputParams = new params();
        inputParams.runTime = Double.parseDouble(args[0]);
        inputParams.startTime = Double.parseDouble(args[1]);
        inputParams.Npatches = Integer.parseInt(args[2]);
        inputParams.S = Integer.parseInt(args[3]);
        inputParams.I = Integer.parseInt(args[4]);
        inputParams.nu = Double.parseDouble(args[5]);
        inputParams.c = Double.parseDouble(args[6]);
        inputParams.U = Double.parseDouble(args[7]);
        test.runSim(inputParams);


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


            double diff = genotype_curr_per_patch.get(o1).size() - genotype_curr_per_patch.get(o2).size();
            if(diff == 0) {

                return 0;
            }
            else if(diff > 0) {
                return -1;
            }
            else{
                return 1;
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
