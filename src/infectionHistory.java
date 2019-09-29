import java.io.*;
import java.lang.reflect.Array;
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
    List<Double> fitness = null;

    public infectionHistory () {

        genotype = new ArrayList<>();
        mutationsFromParent = new ArrayList<>();
        birth = new ArrayList<>();
        death = new ArrayList<>();
        parent = new ArrayList<>();
        prevalence = new ArrayList<>();
        patch = new ArrayList<>();
        fitness = new ArrayList<>();

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
    public void logFitness(Double fitness) {

        this.fitness.add(fitness);

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

    public Double getFitness(int index) {

        return this.fitness.get(index);
    }

    public Integer getParent(Integer genotype) {

        int index = this.genotype.indexOf(genotype);
        return this.parent.get(index);
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
        this.fitness = new ArrayList<>();
        this.fitness.addAll(history.fitness);
    }


    public List<Double> getCompleteBirths() {
        return this.birth;
    }
    public List<Double> getCompleteDeaths() {
        return this.death;
    }


//    public byte[] serializaData() {
//
//        ArrayList<Object> x = new ArrayList<>();
//        x.add(genotype);
//        x.add(mutationsFromParent);
//        x.add(birth);
//        x.add(death);
//        x.add(parent);
//        x.add(prevalence);
//        x.add(patch);
//        x.add(fitness);
//
//        byte[] bytes = null;
//        try {
//            //FileOutputStream fileOut = new FileOutputStream("history.ser");
//            ByteArrayOutputStream bos = new ByteArrayOutputStream();
//            //bos.writeTo(fileOut);
//            ObjectOutputStream out = new ObjectOutputStream(bos);
//            out.flush();
//            out.writeObject(x);
//            bytes = bos.toByteArray();
//            out.close();
//            bos.close();
//            //System.out.println("Serialized data is saved in bytes "+bytes.length);
//            return bytes;
//        } catch (IOException i) {
//            i.printStackTrace();
//        }
//
//        return null;
//
//
//    }
//
//    public void deserializaData(byte[] bytes) {
//
//        //ArrayList<Object> deserialized = null;
//        // Deserialize in to new class object
//        Object deserialized = null;
//
//        try {
//            ByteArrayInputStream bis = new ByteArrayInputStream(bytes);
//
//            //FileInputStream fileIn = new FileInputStream("history.ser");
//
//            ObjectInputStream in = new ObjectInputStream(bis);
//            deserialized = in.readObject();
//            in.close();
//            bis.close();
//        }
//        catch(IOException i){
//            i.printStackTrace();
//        } catch(ClassNotFoundException c){
//            c.printStackTrace();
//        }
//
//        ArrayList<Object> x = (ArrayList<Object>)deserialized;
//        this.genotype = ((ArrayList<Integer>)x.get(0));
//        this.mutationsFromParent = ((ArrayList<Integer>)x.get(1));
//        this.birth = ((ArrayList<Double>)x.get(2));
//        this.death = ((ArrayList<Double>)x.get(3));
//        this.parent = ((ArrayList<Integer>)x.get(4));
//        this.prevalence = ((ArrayList<List<Integer>>)x.get(5));
//        this.patch = ((ArrayList<Integer>)x.get(6));
//        this.fitness = ((ArrayList<Double>)x.get(7));
//
//
//    }
//
//



}
