// #include <string>
// #include <vector>
#include <cstdio>

#include "3.h"
using namespace std;

void readF(const string &configFile, FileName &filename);
void read0(const FileName &filename, Data &data);
void read1(const FileName &filename, parameter1 &p1);
void read3(const FileName &filename, parameter3 &p3);
void read4(const FileName &filename, parameter4 &p4);
void read5(const FileName &filename, parameter5 &p5);
//void read2(const FileName &filename, parameter2 &p2);
void read2(const FileName &filename, parameter1 &p1, parameter2 &p2_WW, parameter2 &p2_BB, parameter2 &p2_WB);
bool ReadDataFrameForResume(const std::string &path,
                            double target_time_fs,
                            double dt_fs,
                            bool use_last_frame,
                            Data &data);
bool ReadEt0ResumeState(const std::string &et0_path,
                        double target_time_fs,
                        double dt_fs,
                        bool use_last_row,
                        const parameter4 &p4,
                        ResumeState &st);
std::string MakeEtGroupPath(const std::string &et_path, int group_index);
void Dimensionalless(Data &data, parameter1 &pr1, parameter2 &pr2_WW, parameter2 &pr2_BB, parameter2 &pr2_WB);
void IntialData(Data &data, parameter1 &p1);
void LJ_potential(Data &data, const parameter1 &pr1, vector<double> &U_atom);
//void BeW_potential(const parameter1 &pr1, const parameter2 &pr2, Data &data, vector<double> &U_atom);
void BeW_potential(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,vector<double> &U_atom);

void lcl0(Data&data,Cell_List&cl,const parameter1 &pr1,
          const parameter2 &pr2_WW,
          const parameter2 &pr2_BB,
          const parameter2 &pr2_WB);
void lcl1(Data&data,Cell_List&cl,const parameter1 &pr1, const parameter2 &pr2,vector<double> &U_atom);
void lcl2(Data&data,Cell_List&cl,const parameter1 &pr1, 
            const parameter2 &pr2_WW,
            const parameter2 &pr2_BB,
            const parameter2 &pr2_WB,
            vector<double> &U_atom);
void Force_BeW_ABOP_LCL(const parameter1 &pr1,
                        const parameter2 &pr2_WW,
                        const parameter2 &pr2_BB,
                        const parameter2 &pr2_WB,
                        Data &data, Cell_List &cl, vector<double> &U_atom);
void Force_Current(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data, Cell_List &cl, vector<double> &U_atom);
bool SetForceModelByName(const std::string &name);
const char *GetForceModelName();
int ProjectileSpawnCount_Current(const parameter5 &p5, const InjectionRuntime &ir, const Data &data);
Matrix31 ProjectileVelocity_Current(const parameter5 &p5, const parameter1 &pr1, const Data &data, const Atom &a_template);
Matrix31 ProjectilePosition_Current(const parameter5 &p5, InjectionRuntime &ir, const Data &data, int local_index);
bool MaybeInjectProjectiles(const parameter5 &p5,
                           const parameter1 &pr1,
                           Data &data,
                           InjectionRuntime &ir,
                           bool verbose = true);
void BeW_potential2(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,vector<double> &U_atom);
void LJ_evolution(const parameter1 &pr1, Data &Data0,vector<double> &U_atom);

/*void BeW_evolution(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,vector<double> &U_atom);
                   */
void BeW_evolution1(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,Cell_List&cl,
                   vector<double> &U_atom);
void BeW_evolution2(const parameter1 &pr1,
                   const parameter2 &pr2_WW,
                   const parameter2 &pr2_BB,
                   const parameter2 &pr2_WB,
                   Data &data,Cell_List&cl,
                   vector<double> &U_atom);
void energy(Data &data, const parameter1 &pr1, vector<double> &U_atom);
void InitEnergyFile(const FileName &filename);
void OutputData(const Data &data, const FileName &filename, const parameter1 &pr1);
void OutputEnergy(const Data &data, const FileName &filename, const parameter1 &pr1); 

class OutputManager
{
public:
    OutputManager();
    ~OutputManager();

    bool Open(const FileName &filename);
    bool Rotate(const FileName &filename);
    void WriteData(const Data &data, const parameter1 &pr1);
    void WriteEnergy(const Data &data, const parameter1 &pr1);
    void Flush();
    void Close();

private:
    FILE *data_fp_;
    std::vector<FILE *> et_group_fp_;
    std::vector<std::string> et_group_name_;
    std::vector<std::string> species_name_cache_;
    std::string et_base_path_;
    bool et_groups_initialized_;
};
