import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jayna on 18/05/2018.
 */
public class infectionHistory implements Serializable {

    List<Integer> genotype = null;
    List<Integer> mutationsFromParent = null;
    List<Double> birth = null;
    List<Double> death = null;
    List<Integer> parent = null;
    List<List<Integer>> prevalence = null;
    List<Integer> patch = null;
    List<Double> beta = null;
    List<Integer> clade = null;


    public infectionHistory () {

        genotype = new ArrayList<>();
        clade = new ArrayList<>();
        mutationsFromParent = new ArrayList<>();
        birth = new ArrayList<>();
        death = new ArrayList<>();
        parent = new ArrayList<>();
        prevalence = new ArrayList<>();
        patch = new ArrayList<>();
        beta = new ArrayList<>();

    }

    public void logClade(Integer clade) {

        this.clade.add(clade);

    }
    public void logGenotype(Integer genotype) {

        this.genotype.add(genotype);

    }


    public void logMutations(Integer mutations) {

        this.mutationsFromParent.add(mutations);

    }

    public void logBirth(Double birth) {

        this.birth.add(birth);
    }

    public void logDeath(Double death) {

        this.death.add(death);
    }

    public void logParent(Integer parent) {

        this.parent.add(parent);
    }

    public void logPatch(Integer patch) {

        this.patch.add(patch);
    }




    public void logPrevalence(Integer prevalence, int index) {

        this.prevalence.get(index).add(prevalence);
    }
    public void logBeta(Double beta) {

        this.beta.add(beta);

    }

    public int getId(int index) {

        return this.genotype.get(index);
    }

    public Double getBirth(int index) {

        return this.birth.get(index);
    }


    public Double getBirth(Integer genotype) {

        int index = this.genotype.indexOf(genotype);
        return this.birth.get(index);
    }

    public int getPrevalenceAtTimet(Integer genotype, int time) {

        int index = this.genotype.indexOf(genotype);
        int prev_at_t = this.prevalence.get(index).get(time);
        return prev_at_t;
    }

    public Integer getMutationsFromParent(Integer genotype) {


        int index = this.genotype.indexOf(genotype);
        return this.mutationsFromParent.get(index);
    }

    public Integer getMutationsFromParent(Integer genotype, Integer parent) {

        if(genotype.equals(parent)) {
            return new Integer(0);
        }
        else {
            Integer mutations = getMutationsFromParent(genotype);
            Integer parent_i = getParent(genotype);


            while (!parent_i.equals(parent)) {

                mutations += getMutationsFromParent(parent_i);
                parent_i = new Integer(getParent(parent_i));

            }

            return mutations;
        }
    }

    public Double getDeath(int index) {

        return this.death.get(index);
    }

    public Integer getPatch(int index) {

        return this.patch.get(index);
    }

    public Integer getClade(int index) {

        return this.clade.get(index);
    }

    public Integer getClade(Integer genotype) {

        int index = this.genotype.indexOf(genotype);
        return this.clade.get(index);
    }

    public Double getBeta(int index) {

        return this.beta.get(index);
    }


    public Double getBeta(Integer genotype) {

        int index = this.genotype.indexOf(genotype);
        return this.beta.get(index);
    }

    public Integer getParent(Integer genotype) {

        int index = this.genotype.indexOf(genotype);
        return this.parent.get(index);
    }

    public void setBeta(Integer genotype, Double beta) {

        int index = this.genotype.indexOf(genotype);
        this.beta.set(index, beta);

    }

    public void copyInfectionHistory(infectionHistory history) {

        this.genotype = new ArrayList<>();
        this.genotype.addAll(history.genotype);
        this.prevalence = new ArrayList<>();
        this.prevalence.addAll(history.prevalence);
        this.death = new ArrayList<>();
        this.death.addAll(history.death);
        this.birth = new ArrayList<>();
        this.birth.addAll(history.birth);
        this.patch = new ArrayList<>();
        this.patch.addAll(history.patch);
        this.mutationsFromParent = new ArrayList<>();
        this.mutationsFromParent.addAll(history.mutationsFromParent);
        this.beta = new ArrayList<>();
        this.beta.addAll(history.beta);
        this.clade = new ArrayList<>();
        this.clade.addAll(history.clade);
    }


    public List<Double> getCompleteBirths() {
        return this.birth;
    }
    public List<Double> getCompleteDeaths() {
        return this.death;
    }





}
