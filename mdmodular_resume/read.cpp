#include <iostream>
// #include <random>
#include <sstream>
#include <string>
// #include <cmath>
// #include <fstream>
#include <vector>
#include <algorithm>
#include <cctype>
// #include <iomanip>

#include "3.h"

//---------------------------------------------------------------------------
// BasicData_filename,Data_filename,Ut_file,Kt_file
void readF(const std::string &configFile, FileName &fn)
{
    fn.parameter3_filename = "parameter.p3";
    fn.parameter4_filename = "parameter.p4";
    fn.parameter5_filename = "parameter.p5";

    std::ifstream fin(configFile);
    if (!fin)
    {
        std::cerr << "can't open file: " << configFile << std::endl;
        return;
    }

    string key, x, line;

    while (getline(fin, line))
    {
        if (line.empty() ||
            line[0] == '#' ||
            (line.size() > 1 && line[0] == '/' && line[1] == '/'))
            continue;

        std::stringstream ss(line);
        ss >> key >> x;
        // parameter                           parameter.p1
        // BCC

        // First_Data_file                          /BasicData/BCC.md3
        // Data_file                                /Data/BCC.md3
        // Et_file                             /Data/BCC.Et.md3
        if (key == "")
            continue; // skip space
        if (x == "")
            continue; // 

        if (key == "parameter")
        {
            fn.parameter_filename = x;
        }
        else if (key == "parameter2")
        {
            fn.parameter2_filename = x;
        }
        else if (key == "parameter3")
        {
            fn.parameter3_filename = x;
        }
        else if (key == "parameter4")
        {
            fn.parameter4_filename = x;
        }
        else if (key == "parameter5")
        {
            fn.parameter5_filename = x;
        }
        else if (key == "First_Data_file")
        {
            fn.BasicData_filename = x;
        }
        else if (key == "Data_file")
        {
            fn.Data_filename = x;
        }
        else if (key == "Et_file")
        {
            fn.Et_file = x;
        }
    }
}

void read0(const FileName &filename, Data &data)
{
    ifstream fin(filename.BasicData_filename);
    if (!fin)
    {
        cerr << "Can't open file: " << filename.BasicData_filename << endl;
        return;
    }

    // 1. read n
    fin >> data.n;
    string dummy;
    getline(fin, dummy);

    // 2. read time & energy
    string tmp;
    fin >> tmp >> data.t >> tmp >> tmp >> data.E;
    getline(fin, dummy);

    // 3. read Box
    fin >> tmp;
    fin >> data.Box.a00 >> data.Box.a01 >> data.Box.a02 >> data.Box.a10 >> data.Box.a11 >> data.Box.a12 >> data.Box.a20 >> data.Box.a21 >> data.Box.a22;
    getline(fin, dummy);

    // 4. read atoms
    data.atoms.resize(data.n);

    for (int i = 0; i < data.n; ++i)
    {
        std::string atom_line;
        if (!std::getline(fin, atom_line))
            break;
        if (atom_line.empty())
        {
            --i;
            continue;
        }

        std::stringstream ss(atom_line);
        Atom &a = data.atoms[i];
        a.identity = 0;
        ss >> a.name >> a.r.a00 >> a.r.a10 >> a.r.a20
           >> a.p.a00 >> a.p.a10 >> a.p.a20
           >> a.f.a00 >> a.f.a10 >> a.f.a20;
        if (!ss)
        {
            std::cerr << "Bad atom line in " << filename.BasicData_filename << " at atom " << i << "\n";
            break;
        }
        if (!(ss >> a.identity))
            a.identity = 0; // backward-compatible: old files have no identity column
        a.atom_type = AtomTypeFromName(a.name);
    }
}

void read1(const FileName &filename, parameter1 &p1)
{
    ifstream fin(filename.parameter_filename);
    if (!fin)
    {
        std::cerr << "can't open file: " << filename.parameter_filename << std::endl;
        return;
    }

    std::string line;
    std::string key;
    char eq = 0;
    double value = 0.0;

    // 
    double endtime_fs = 0.0;
    double steps_space_fs = 0.0;

    while (std::getline(fin, line))
    {
        // delete comment#
        std::size_t pos = line.find('#');
        if (pos != std::string::npos)
            line = line.substr(0, pos);

        // delete space
        auto not_space = [](int ch)
        { return !std::isspace(ch); };

        if (!line.empty())
        {
            line.erase(line.begin(), std::find_if(line.begin(), line.end(), not_space));
            line.erase(std::find_if(line.rbegin(), line.rend(), not_space).base(), line.end());
        }

        if (line.empty())
            continue;

        std::stringstream ss(line);
        ss >> key >> eq >> value;
        if (!ss || eq != '=')
            continue;

        if (key == "dt")
        {
            p1.dt = value; // fs
        }
        else if (key == "endtime")
        {
            endtime_fs = value; // fs
        }
        else if (key == "steps_space_fs")
        {
            // space time（fs）
            steps_space_fs = value;
        }
        else if (key == "T_end")
        {
            p1.T_end = value;
        }
        else if (key == "dT")
        {
            p1.dT = value;
        }
        else if (key == "T_init")
        {
            p1.TT = value;
        }
        else if (key == "T_time")
        {
            p1.T_time = value;
        }
        else if (key == "output_Data_time")
        {
            p1.output_Data_time = value;
        }
        else if (key == "output_Et_time")
        {
            p1.output_Et_time = value;
        }
        else if (key == "mass_W")
        {
            p1.mw = value;
        }
        else if (key == "mass_Be")
        {
            p1.mb = value;
        }
        else if (key == "kb")
        {
            p1.kb = value;
        }
            else if (key == "E0")
        {
            p1.E0 = value;
        }
        else if (key == "epsilon")
        {
            p1.epsilon = value;
        }
        else if (key == "sigma")
        {
            p1.sigma = value;
        } 
        else if (key == "s0")
        {
            p1.s0 = value;
        } 
        else if (key == "ps0")
        {
            p1.ps0 = value;
        }
        else if (key == "Q")
        {
            p1.Q = value;
        }
        else if (key == "xi")
        {
            p1.xi = value;
        }
        else if (key == "H0")
        {
            p1.H0 = value;
        }

    }

    // culculate steps and transpose to int
    if (p1.dt > 0.0 && endtime_fs > 0.0)
    {
        p1.steps = static_cast<int>(std::round(endtime_fs / p1.dt));
    }

    // steps_space to int
    if (p1.dt > 0.0 && steps_space_fs > 0.0)
    {
        p1.steps_space = static_cast<int>(std::round(steps_space_fs / p1.dt));
    }

    fin.close();
}

void read3(const FileName &filename, parameter3 &p3)
{
    std::string path = filename.parameter3_filename.empty() ? "parameter.p3" : filename.parameter3_filename;
    std::ifstream fin(path);
    if (!fin)
    {
        std::cerr << "can't open file: " << path << std::endl;
        return;
    }

    std::string line;
    while (std::getline(fin, line))
    {
        std::size_t pos = line.find('#');
        if (pos != std::string::npos)
            line = line.substr(0, pos);

        auto not_space = [](int ch)
        { return !std::isspace(ch); };
        if (!line.empty())
        {
            line.erase(line.begin(), std::find_if(line.begin(), line.end(), not_space));
            line.erase(std::find_if(line.rbegin(), line.rend(), not_space).base(), line.end());
        }
        if (line.empty())
            continue;

        std::string key;
        char eq = 0;
        std::string value_str;
        std::stringstream ss(line);
        ss >> key >> eq >> value_str;
        if (!ss || eq != '=')
            continue;

        if (key == "extra_steps_max")
            p3.extra_steps_max = std::stoi(value_str);
        else if (key == "plateau_blocks")
            p3.plateau_blocks = std::stoi(value_str);
        else if (key == "output_rotate_every_temp_up")
            p3.output_rotate_every_temp_up = std::stoi(value_str);
        else if (key == "output_rotate_every_temp_down")
            p3.output_rotate_every_temp_down = std::stoi(value_str);
        else if (key == "output_flush_every_writes")
            p3.output_flush_every_writes = std::stoi(value_str);
        else if (key == "verlet_skin")
            p3.verlet_skin = std::stod(value_str);
    }
}

void read4(const FileName &filename, parameter4 &p4)
{
    std::string path = filename.parameter4_filename.empty() ? "parameter.p4" : filename.parameter4_filename;
    std::ifstream fin(path);
    if (!fin)
    {
        std::cerr << "can't open file: " << path << std::endl;
        return;
    }

    std::string line;
    while (std::getline(fin, line))
    {
        std::size_t pos = line.find('#');
        if (pos != std::string::npos)
            line = line.substr(0, pos);

        auto not_space = [](int ch)
        { return !std::isspace(ch); };
        if (!line.empty())
        {
            line.erase(line.begin(), std::find_if(line.begin(), line.end(), not_space));
            line.erase(std::find_if(line.rbegin(), line.rend(), not_space).base(), line.end());
        }
        if (line.empty())
            continue;

        std::string key;
        char eq = 0;
        std::string value_str;
        std::stringstream ss(line);
        ss >> key >> eq >> value_str;
        if (!ss || eq != '=')
            continue;

        if (key == "continue_run")
            p4.continue_run = std::stoi(value_str);
        else if (key == "resume_use_target_time")
            p4.resume_use_target_time = std::stoi(value_str);
        else if (key == "resume_time_fs")
            p4.resume_time_fs = std::stod(value_str);
        else if (key == "resume_require_exact_time")
            p4.resume_require_exact_time = std::stoi(value_str);
        else if (key == "continue_write_next_index")
            p4.continue_write_next_index = std::stoi(value_str);
        else if (key == "force_model")
            p4.force_model = value_str;
    }
}

void read5(const FileName &filename, parameter5 &p5)
{
    std::string path = filename.parameter5_filename.empty() ? "parameter.p5" : filename.parameter5_filename;
    std::ifstream fin(path);
    if (!fin)
    {
        std::cerr << "can't open file: " << path << std::endl;
        return;
    }

    std::string line;
    while (std::getline(fin, line))
    {
        std::size_t pos = line.find('#');
        if (pos != std::string::npos)
            line = line.substr(0, pos);

        auto not_space = [](int ch)
        { return !std::isspace(ch); };
        if (!line.empty())
        {
            line.erase(line.begin(), std::find_if(line.begin(), line.end(), not_space));
            line.erase(std::find_if(line.rbegin(), line.rend(), not_space).base(), line.end());
        }
        if (line.empty())
            continue;

        std::string key;
        char eq = 0;
        std::string value;
        std::stringstream ss(line);
        ss >> key >> eq >> value;
        if (!ss || eq != '=')
            continue;

        if (key == "enable_injection")
            p5.enable_injection = std::stoi(value);
        else if (key == "incident_species")
            p5.incident_species = value;
        else if (key == "incident_identity")
            p5.incident_identity = std::stoi(value);
        else if (key == "inject_start_time_fs")
            p5.inject_start_time_fs = std::stod(value);
        else if (key == "inject_interval_fs")
            p5.inject_interval_fs = std::stod(value);
        else if (key == "inject_max_particles")
            p5.inject_max_particles = std::stoi(value);
        else if (key == "inject_batch_size")
            p5.inject_batch_size = std::stoi(value);
        else if (key == "inject_z_margin")
            p5.inject_z_margin = std::stod(value);
        else if (key == "inject_x_margin")
            p5.inject_x_margin = std::stod(value);
        else if (key == "inject_y_margin")
            p5.inject_y_margin = std::stod(value);
        else if (key == "inject_speed")
            p5.inject_speed = std::stod(value);
        else if (key == "inject_use_fixed_seed")
            p5.inject_use_fixed_seed = std::stoi(value);
        else if (key == "inject_seed")
            p5.inject_seed = static_cast<unsigned int>(std::stoul(value));
    }
}

void read2(const FileName &filename, parameter1 &p1, parameter2 &pr2_WW, parameter2 &pr2_BB, parameter2 &pr2_WB)
{
    ifstream fin(filename.parameter2_filename);
    if (!fin)
    {
        std::cerr << "can't open file: " << filename.parameter2_filename << std::endl;
        return;
    }

    std::string line;

    auto reset_abop = [](parameter2 &p)
    {
        p.D0 = 0.0;
        p.r0 = 0.0;
        p.beta = 0.0;
        p.S = 0.0;
        p.gamma = 0.0;
        p.c = 0.0;
        p.d = 0.0;
        p.h = 0.0;
        p.R = 0.0;
        p.D = 0.0;
        p.mu = 0.0;
        p.rf = 0.0;
        p.bf = 0;
    };

    reset_abop(pr2_WW);
    reset_abop(pr2_BB);
    reset_abop(pr2_WB);

    auto lower_copy = [](std::string s)
    {
        for (size_t i = 0; i < s.size(); ++i)
            s[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(s[i])));
        return s;
    };

    while (std::getline(fin, line))
    {
        // delete comment#
        std::size_t pos = line.find('#');
        if (pos != std::string::npos)
            line = line.substr(0, pos);

        // delete space
        auto not_space = [](int ch)
        { return !std::isspace(ch); };

        if (!line.empty())
        {
            line.erase(line.begin(), std::find_if(line.begin(), line.end(), not_space));
            line.erase(std::find_if(line.rbegin(), line.rend(), not_space).base(), line.end());
        }

        if (line.empty())
            continue;

        // such as ABB.D0 = xx
        std::string key;
        char eq = 0;
        double value = 0.0;

        std::stringstream ss(line);
        ss >> key;

        if (!ss)
            continue;
        const std::string key_l = lower_copy(key);

        // Legacy mass keys without dot, placed in .pr2 file for convenience.
        if (key_l == "mass_w" || key_l == "massw")
        {
            ss >> eq >> value;
            if (ss && eq == '=')
                p1.mw = value;
            continue;
        }
        if (key_l == "mass_be" || key_l == "massbe")
        {
            ss >> eq >> value;
            if (ss && eq == '=')
                p1.mb = value;
            continue;
        }

        // 拆分出前缀和参数名
        std::size_t dot_pos = key.find('.');
        if (dot_pos == std::string::npos)
        {
            //  if not PREFIX.name  then skip
            continue;
        }

        std::string prefix = key.substr(0, dot_pos);       // ABB / AWW / ABW
        std::string pname = key.substr(dot_pos + 1);       // D0 / r0 / ...

        // Allow masses to live in the potential parameter file, e.g.
        // mass.W = 183.84
        // mass.Be = 9.012
        // Also accepts mass_w / mass_be style keys in this file.
        const std::string prefix_l = lower_copy(prefix);
        const std::string pname_l = lower_copy(pname);
        if (prefix_l == "mass")
        {
            ss >> eq >> value;
            if (!ss || eq != '=')
                continue;

            if (pname_l == "w" || pname_l == "tungsten")
                p1.mw = value;
            else if (pname_l == "be" || pname_l == "beryllium")
                p1.mb = value;
            continue;
        }

        ss >> eq >> value;
        if (!ss || eq != '=')
            continue;

        parameter2 *target = nullptr;
        if (prefix == "AWW" || prefix == "Aww" || prefix == "aww")
        {
            target = &pr2_WW; // W-W
        }
        else if (prefix == "ABB" || prefix == "Abb" || prefix == "abb")
        {
            target = &pr2_BB; // Be-Be
        }
        else if (prefix == "AWB" || prefix == "Awb" || prefix == "awb")
        {
            target = &pr2_WB; // Be-W
        }
        else
        {
            continue;
        }

        if (pname == "D0")
            target->D0 = value;
        else if (pname == "r0")
            target->r0 = value;
        else if (pname == "beta")
            target->beta = value;
        else if (pname == "S")
            target->S = value;
        else if (pname == "gamma")
            target->gamma = value;
        else if (pname == "c")
            target->c = value;
        else if (pname == "d")
            target->d = value;
        else if (pname == "h")
            target->h = value;
        else if (pname == "R")
            target->R = value;
        else if (pname == "D")
            target->D = value;
        else if (pname == "2mu")
            //target->mu = 0.5 * value; // convert 2mu -> mu
            target->mu = value; // convert 2mu -> mu
        //else if (pname == "mu")
            //target->mu = value;
        else if (pname == "rf")
            target->rf = value;
        else if (pname == "bf")
            target->bf = value;
    }

    fin.close();

    std::cout << "[read2] masses in use: W=" << p1.mw << " Be=" << p1.mb
              << " (parameter2 file may override parameter1)\n";
}
