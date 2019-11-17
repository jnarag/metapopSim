/**
 * Created by jayna on 17/05/2018.
 */

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.jet.random.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

import org.eclipse.collections.api.list.primitive.MutableIntList;
import org.eclipse.collections.impl.factory.Iterables;
import org.eclipse.collections.impl.factory.primitive.IntLists;
import sun.util.resources.cldr.ca.CurrencyNames_ca;
import utils.*;


public class patchSim {

    params params;
    int genotype;
    infectionHistory patch_history;
    Map<Integer, List<Integer>> genotype_curr_per_patch;
    List<Integer> occupyPatchList;
    DoubleSummaryStatistics patchLifespan;


    int total_infected;
    int total_genotypes;
    double threshold;
    double r_i;

    public patchSim() {

        patch_history = new infectionHistory();
        genotype = 1;
        genotype_curr_per_patch = new HashMap<>();
        occupyPatchList = new ArrayList<>();
        patchLifespan = new DoubleSummaryStatistics();

        total_genotypes = 0;
        total_infected = 0;
        r_i = params.r;

    }


    public void runSim(params params) {


        double tau = params.tau;
        boolean writeOutput = false;

        String outputfile = "patchSim_extinction_" + params.nu + "_col_" + params.c +
                "_beta_" + params.beta + "_Npatches_" + params.Npatches +
                "_patchSize_" + (params.S + params.I + ".txt");


        String summaryFile = "summary_patchSim_extinction_" + params.nu + "_col_" + params.c +
                "_beta_" + params.beta + "_Npatches_" + params.Npatches +
                "_patchSize_" + (params.S + params.I + ".txt");

        FileWriter writer1 = null;
        FileWriter writer2 = null;
        if(writeOutput) {
            try {
                writer1 = new FileWriter(new File(outputfile));
                writer1.write("Time\tPatch_no\tTotal_infected\tUnique_genotypes\n");

                writer2 = new FileWriter(new File(summaryFile));
                writer2.write("Time\tTotal_infected\tUnique_genotypes\tOccupied_patches\tGenealogical_diversity\tTMRCA\n");
            } catch (IOException e) {
                e.printStackTrace();
            }

        }

        double t_curr = tau;
        double t_max = params.runTime;
        int N = params.S + params.I;

        DoubleArrayList betweenPatchRates = new DoubleArrayList();
        DoubleArrayList withinPatchRates = new DoubleArrayList();


        double[] timeRefactory_curr_i = new double[params.Npatches];
        threshold = 0.75*N;


        //initializing patches with number of infected and susceptible cells

        Arrays.fill(timeRefactory_curr_i, -5);


        double colonize = params.c;
        double extinct = params.nu;
        String whModel = params.whModel;


        int init_Npatches = 1; //(int)Math.floor(params.Npatches*(colonize/extinct));

        System.out.println("Number of occupied patches at equilibrium "+ (int)Math.floor(params.Npatches*(colonize/extinct)));

        System.out.println("Number of cells in the metapopulation "+ (Math.log10(1.*params.Npatches*N)));

        initialiseHistory(patch_history, params.I, init_Npatches, t_max);


        for(Integer genotype: patch_history.genotype) {
            int genotype_index = patch_history.genotype.indexOf(genotype);

            int patch = patch_history.getPatch(genotype_index);
            occupyPatchList.add(patch);

            int n_infected = patch_history.prevalence.get(genotype_index).get(0);
            List<Integer> Y_temp = new ArrayList<>(Arrays.asList(new Integer[n_infected]));
            Collections.fill(Y_temp, genotype);

            genotype_curr_per_patch.put(patch, Y_temp);

            total_infected += n_infected;
            total_genotypes += 1;

        }


        String sim = "ext_"+ extinct+"_col_"+colonize+"_Npatch_"+ params.Npatches +"_PatchSize_"+(N);

        int N_occupied = occupyPatchList.size();


        for (int t = 1; t < (t_max / tau); t++) {


//            //determine what is present currently across all patches
            update_curr_in_body(t);


            if (t_curr % 20 == 0) {

                updatePatchStatistics(writer1, t_curr, writeOutput);

                List<Double> globalDiversity = updateDiversity(false);
                double tmrca = globalDiversity.get(1);

                System.out.println("t_curr \t" + t_curr + "\ttotal_infected: \t" + Math.log10(total_infected) + "\ttotal_genotypes: \t" + total_genotypes +
                        "\tpatch_occupied: \t" + N_occupied + "\t" + (tmrca));


                if (writeOutput) {


                    writeOutput(writer2, t_curr, total_infected, total_genotypes, N_occupied, globalDiversity);
                }

            }


            betweenPatchRates.add(extinct * N_occupied);
            betweenPatchRates.add(colonize * (params.Npatches - N_occupied));

            int noOfRates = betweenPatchRates.size();

            Poisson poisson;

            int N_occupied_copy = N_occupied;
            for(int event = 0; event < noOfRates; event++) {

                int j;

                poisson = new Poisson(betweenPatchRates.get(event)*tau, params.randomGenerator);
                int num = poisson.nextInt();

                if(num == 0) continue;

                //extinction

                switch(event) {
                    case 0:

                        int min = Math.min(N_occupied, num);

                        j = 0;

                        while (j < min) {

                            // patch extinction

                            int patch_index = Uniform.staticNextIntFromTo(0, (N_occupied - 1));
                            Integer patch = occupyPatchList.get(patch_index);

                            // track infection history (genotype_curr should be set to 0)

                            Integer g = genotype_curr_per_patch.get(patch).get(0);


                            int index = patch_history.genotype.indexOf(g);

                            //patch_history.prevalence.get(index).set(t, 0);
                            patch_history.death.set(index, t_curr);

                            total_infected -= genotype_curr_per_patch.get(patch).size();
                            total_genotypes -= 1;


                            genotype_curr_per_patch.remove(patch);
                            occupyPatchList.remove(patch_index);



                            if(occupyPatchList.size() != genotype_curr_per_patch.size()) {
                                System.out.println("trouble");
                            }
                            //thresholdList.remove((Integer) patch_index);
                            timeRefactory_curr_i[patch_index] = t_curr;
                            N_occupied -= 1;

                            double patchBirth = patch_history.birth.get(index);

                            patchLifespan.accept(t_curr-patchBirth);



                            j++;

                        }
                        break;

                    // patch colonisation


                    case 1:

                        min = Math.min((params.Npatches-N_occupied), num);

                        j = 0;

                        while (j < min) {

                            j++;
                            int new_patch = Uniform.staticNextIntFromTo(0, params.Npatches - 1);

                            if(occupyPatchList.contains(new_patch)) {
                                while(occupyPatchList.contains(new_patch)) {
                                    new_patch = Uniform.staticNextIntFromTo(0, params.Npatches - 1);
                                }
                            }

                            if (total_infected == 0) {
                                break;
                            }

                            if (t_curr - timeRefactory_curr_i[new_patch] < 5) {
                                continue;
                            }

                            int s = Uniform.staticNextIntFromTo(0, N_occupied - 1);

                            Integer source = occupyPatchList.get(s);


                            int newInfections = Uniform.staticNextIntFromTo(1, 3);


                            Integer colonizing_genotype = genotype_curr_per_patch.get(source).get(0);  //curr_in_body.get(r);

                            List<Integer> new_Y = new ArrayList<>(Arrays.asList(new Integer[newInfections]));
                            Collections.fill(new_Y, genotype);
                            genotype_curr_per_patch.put(new_patch, new_Y);

                            occupyPatchList.add(new_patch);


                            total_infected += newInfections;
                            total_genotypes += 1;

                            patch_history.logBirth(t_curr);
                            patch_history.logDeath(Double.POSITIVE_INFINITY);
                            patch_history.logGenotype(genotype);
                            patch_history.logParent(colonizing_genotype);
                            patch_history.logPatch(new_patch);

                            List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau)]));
                            Collections.fill(new_prevalence, 0);
                            new_prevalence.set(t, newInfections);
                            patch_history.prevalence.add(new_prevalence);
                            genotype++;


                            if(genotype_curr_per_patch.size()!= occupyPatchList.size()) {
                                System.out.println("trouble!");
                            }
                            N_occupied += 1;

                        }
                        break;
                }

            }

            betweenPatchRates.clear();


            for (int p = 0; p < N_occupied; p++) {

                Integer patch = occupyPatchList.get(p);
                List<Integer> genotype_curr = new ArrayList<>();
                genotype_curr.addAll(genotype_curr_per_patch.get(patch));

                int Y_temp = genotype_curr.size();
                int Y_temp_copy = Y_temp;

                if(Y_temp == 0) continue;

                double infection_rate;
                double death;
                double death_rate;

                switch(whModel) {

                    case "log":
                        infection_rate = ((r_i*Y_temp)*(1.0-(Y_temp*1.0/N*1.0)));
                        withinPatchRates.add(infection_rate);
                        break;
                    case "BSI_death":

                        //BSI - beta is scaled by total number of patches. Should this be scaled by number of occupied patches?
                        infection_rate = (params.beta/N_occupied)*(N-Y_temp)*(Y_temp);
                        death = params.death/N_occupied;
                        death_rate = death*Y_temp;
                        withinPatchRates.add(infection_rate);
                        withinPatchRates.add(death_rate);

                        break;
                    case "BSI":

                        infection_rate = (params.beta/N_occupied)*(N-Y_temp)*(Y_temp);
                        withinPatchRates.add(infection_rate);
                        break;


                }

                noOfRates = withinPatchRates.size();

                for(int event = 0; event < noOfRates; event++) {

                    poisson = new Poisson(withinPatchRates.get(event) * tau, params.randomGenerator);
                    int num = poisson.nextInt();

                    if(num == 0) continue;
                    // max number of new infected cells is bounded by number of uninfected cells

                    if(Y_temp == 0) continue;

                    Integer chosen_genotype = genotype_curr.get(0);
                    int genotype_index = patch_history.genotype.indexOf(chosen_genotype);

                    switch (event) {

                        case 0:
                            int min = Math.min((N - Y_temp_copy), num);

                            if(min == 0) continue;

                            int prevalence_curr = patch_history.prevalence.get(genotype_index).get(t - 1);
                            int prevalence_next = patch_history.prevalence.get(genotype_index).get(t);

                            if (prevalence_next > 0) {

                                prevalence_curr = prevalence_next;
                            }

                            List<Integer> temp_list = new ArrayList<>(Arrays.asList(new Integer[min]));
                            Collections.fill(temp_list, chosen_genotype);
                            genotype_curr.addAll(temp_list);

                            patch_history.prevalence.get(genotype_index).set(t, prevalence_curr + min);
                            total_infected += min;

                            genotype_curr_per_patch.replace(patch, genotype_curr);

                            Y_temp += min;

                            break;

                        case 1:

                            min = Math.min(Y_temp_copy, num);

                            if(min == 0) continue;

                            prevalence_curr = patch_history.prevalence.get(genotype_index).get(t - 1);
                            prevalence_next = patch_history.prevalence.get(genotype_index).get(t);


                            if (prevalence_next > 0) {

                                prevalence_curr = prevalence_next;
                            }

                            IntStream.range(0, min).forEach(i -> genotype_curr.remove(chosen_genotype));


                            patch_history.prevalence.get(genotype_index).set(t, prevalence_curr - min);
                            total_infected -= min;

                            genotype_curr_per_patch.replace(patch, genotype_curr);

                            if(genotype_curr.isEmpty()) {
                                occupyPatchList.remove(p);
                                genotype_curr_per_patch.remove(patch);
                                if(genotype_curr_per_patch.size()!= occupyPatchList.size()) {
                                    System.out.println("trouble!");
                                }
                                N_occupied -= 1;
                                timeRefactory_curr_i[patch] = t_curr;
                                total_genotypes -= 1;
                                patch_history.death.set(genotype_index, t_curr);
                                double patchBirth = patch_history.birth.get(genotype_index);

                                patchLifespan.accept(t_curr-patchBirth);

                            }

                            Y_temp -= min;

                            break;
                    }
                }
                withinPatchRates.clear();
            }

            t_curr += tau;

        }

        System.out.println();
        System.out.println("mean: "+patchLifespan.getAverage()+", min: "+patchLifespan.getMin()+ ", max: "+ patchLifespan.getMax());


        if(writeOutput) {
            try {
                writer1.close();
                writer2.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }


    private void initialiseHistory(infectionHistory history, int Y_init, int n_patches, double t_max) {

        int i = 0;
        while(i < n_patches) {

//            if(Y_curr_i[i] > 0) {
            history.logGenotype(genotype);
            history.logBirth(0.0);
            history.logDeath(Double.POSITIVE_INFINITY);
            history.logParent((int) Double.NEGATIVE_INFINITY);
            List<Integer> prevalence_through_time = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / params.tau)]));
            Collections.fill(prevalence_through_time, 0);
            prevalence_through_time.set(0, Y_init);
            history.prevalence.add(prevalence_through_time);
            history.logPatch(i);
            genotype++;
//            }
            i++;
        }
    }

    public infectionHistory getInfectionHistory() {

        return this.patch_history;
    }


    private void writeOutput(FileWriter writer, double t, int total_infected, int genotypes, int n_patches, List<Double> diversity) {

        try {
            writer.write(t+"\t"+total_infected+"\t"+genotypes+"\t"+n_patches+"\t"+diversity.get(0)+"\t"+diversity.get(1) +"\n");
            writer.flush();
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

        double distA = patch_history.getMutationsFromParent(i1, ancestor);
        double distB = patch_history.getMutationsFromParent(i2, ancestor);
        return distA + distB;

    }

    public void updatePatchStatistics(FileWriter writer, double t_curr, boolean writeOutput) {

//        total_infected = 0;
//        total_genotypes = 0;

        for (Map.Entry<Integer, List<Integer>> entry : genotype_curr_per_patch.entrySet()) {

//            for (Integer patch: genotype_curr_per_patch.keySet()) {
            Integer patch = entry.getKey();
            List<Integer> genotype_curr = entry.getValue();

            int Y_i = genotype_curr.size();

            int unique_genotypes = new HashSet<>(genotype_curr).size();
//            total_infected += Y_i;
//            total_genotypes += unique_genotypes;
            int write = -1;

            if (writeOutput == false) {
                continue;
            } else {
                write = 0;
            }

            int check = Y_i + write;


            if (check > 0) {

                try {
                    writer.write(t_curr + "\t" + (patch + 1) + "\t" + Y_i +
                            "\t" + unique_genotypes + "\n");
                    writer.flush();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

        }

    }

    public void update_curr_in_body(int t) {

        for (Map.Entry<Integer, List<Integer>> entry: genotype_curr_per_patch.entrySet()) {
            List<Integer> genotype_curr = entry.getValue();
            try {
                int index = patch_history.genotype.indexOf(genotype_curr.get(0));
                //List<Integer> prevalence = patch_history.prevalence.get(i);
                int prev_curr = patch_history.prevalence.get(index).get(t - 1);
                patch_history.prevalence.get(index).set(t, prev_curr);
            }
            catch (IndexOutOfBoundsException e) {
                e.printStackTrace();
                System.out.println();
            }

        }

    }

    public List<Double> updateDiversity(boolean getDivergence) {

        double diversity = 0.0;
        double tmrca1 = 0.0;
        int sampleCount1 = 0;

        //double geneticDiversity = 0.0;
        double divergence = 0.0;



        List<Integer> weights = new ArrayList<>();


        for(Map.Entry<Integer, List<Integer>> entry: genotype_curr_per_patch.entrySet()) {

            int weight = entry.getValue().size();
            weights.add(weight);

        }

        List<Integer> chosenPatches = new ArrayList<>();
        if(occupyPatchList.size() < 100) {
            chosenPatches.addAll(occupyPatchList);
        }
        else{

            chosenPatches = utils.multinomSamp(weights, occupyPatchList, 100);
        }


        List<Integer> sampled_genotypes = new ArrayList<>();
        chosenPatches.forEach(i -> {
            sampled_genotypes.add(genotype_curr_per_patch.get(i).get(0));

            ;});



        List<Double> results = new ArrayList<>();

        if(sampled_genotypes.size()==0) {
            results.add(diversity);
            results.add(tmrca1);
            //results.add(geneticDiversity);
            if(getDivergence) {
                results.add(divergence);
            }

            return results;
        }


        for (int i = 0; i < sampled_genotypes.size(); i++) {

            Integer vA = sampled_genotypes.get(Uniform.staticNextIntFromTo(0, sampled_genotypes.size()-1));
            Integer vB = sampled_genotypes.get(Uniform.staticNextIntFromTo(0, sampled_genotypes.size()-1));

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
            diversity /= sampleCount1;
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
        inputParams.beta = Double.parseDouble(args[7]);

        test.runSim(inputParams);

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



}
