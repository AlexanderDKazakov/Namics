#include "mesodyn.h"


/* Mesoscale dynamics module written by Daniel Emmery as part of a master's thesis, 2018 */
/* Most of the physics in this module is based on the work of Fraaije et al. in the 1990s  */

Mesodyn::Mesodyn(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Solve_scf*> New_, string name_)
    : Lattice_Access(Lat_[0]),
      name{name_},
      In{In_},
      Lat{Lat_},
      Mol{Mol_},
      Seg{Seg_},
      Sta{Sta_},
      Rea{Rea_},
      Sys{Sys_},
      New{New_},
      D{0.01},
      dt{0.1},
      mean{0},
      stddev{2 * D * sqrt(dt)},
      seed{1},
      seed_specified{false},
      timesteps{100},
      timebetweensaves{1},
      initialization_mode{INIT_HOMOGENEOUS},
      component_no{Sys[0]->SysMolMonList.size()},
      RC{0},
      cn_ratio{0.5},
      component(0),
      flux(0),
      writes{0},
      write_vtk{false},
      order_parameter{0}
  {
  KEYS.push_back("timesteps");
  KEYS.push_back("timebetweensaves");
  KEYS.push_back("delta_t");
  KEYS.push_back("read_pro");
  KEYS.push_back("read_vtk");
  KEYS.push_back("diffusionconstant");
  KEYS.push_back("seed");
  KEYS.push_back("mean");
  KEYS.push_back("stddev");
  KEYS.push_back("cn_ratio");

  // TODO: implement this properly
  // If the user has asked for vtk output
  if (std::find(In[0]->OutputList.begin(), In[0]->OutputList.end(), "vtk") != In[0]->OutputList.end())
    write_vtk = true;

  //Apparently this step takes a ton of time.
  Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[0], writes, timesteps / timebetweensaves));

  set_update_lists();

  // to get the correct KSAM and volume.
  Sys[0]->PrepareForCalculations();

  rho.resize(component_no*M);

  // The right hand side of the minus sign calculates the volume of the boundaries. I know, it's hideous, but it works for all dimensions.
  boundaryless_volume = Sys[0]->volume - ((2 * dimensions - 4) * Lat[0]->MX * Lat[0]->MY + 2 * Lat[0]->MX * Lat[0]->MZ + 2 * Lat[0]->MY * Lat[0]->MZ + (-2 + 2 * dimensions) * (Lat[0]->MX + Lat[0]->MY + Lat[0]->MZ) + pow(2, dimensions));
}

Mesodyn::~Mesodyn() {
    // We only use smart pointers here, they'll take care of deleting themselves when needed.
}

bool Mesodyn::CheckInput(int start) {
  bool success = true;

  /* When adding new items here, please add them to the function prepareOutputFile()*/
  success = In[0]->CheckParameters("mesodyn", name, start, KEYS, PARAMETERS, VALUES);

  if (success) {
    vector<string> options;
    if (GetValue("timesteps").size() > 0) {
      success = In[0]->Get_int(GetValue("timesteps"), timesteps, 1, 100000000, "The number of timesteps should be between 1 and 10000");
    }

    if (GetValue("timebetweensaves").size() > 0) {
      timebetweensaves = In[0]->Get_int(GetValue("timebetweensaves"), timebetweensaves);
    }

    if (GetValue("delta_t").size() > 0) {
      dt = In[0]->Get_Real(GetValue("delta_t"), dt);
    }

    if (GetValue("read_pro").size() > 0) {
      read_filename = In[0]->Get_string(GetValue("read_pro"), read_filename);
      initialization_mode = INIT_FROMPRO;
    }

    if (GetValue("read_vtk").size() > 0) {
      read_filename = In[0]->Get_string(GetValue("read_vtk"), read_filename);
      if (read_filename.find(".vtk") != string::npos) {
        cerr << "Mesodyn will add the component number and extension by itself (in that order), please format the remainder of the filename accordingly." << endl;
        exit(0);
      }
      initialization_mode = INIT_FROMVTK;
    }

    if (GetValue("diffusionconstant").size() > 0) {
      D = In[0]->Get_Real(GetValue("diffusionconstant"), D);
    }

    if (GetValue("seed").size() > 0) {
      seed_specified = true;
      seed = In[0]->Get_Real(GetValue("seed"), seed);
    }

    if (GetValue("mean").size() > 0) {
      mean = In[0]->Get_Real(GetValue("mean"), mean);
    }

    if (GetValue("stddev").size() > 0) {
      stddev = In[0]->Get_Real(GetValue("stddev"), stddev);
    }

    if (GetValue("cn_ratio").size() > 0) {
      cn_ratio = In[0]->Get_Real(GetValue("cn_ratio"), cn_ratio);
    }
  }

  return success;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {

  //   Initialize densities
  initial_conditions();

  //Prepare IO
  set_filename();

  cout << "Mesodyn is all set, starting calculations.." << endl;

  // Do one explicit step before starting the crank nicolson scheme
  explicit_start();

  // Prepare callback functions for SolveMesodyn in Newton
  function<Real*()> solver_callback = bind(&Mesodyn::solve_crank_nicolson, this);
  function<void(Real*,size_t)> loader_callback = bind(&Mesodyn::load_alpha, this, std::placeholders::_1, std::placeholders::_2);

  /**** Main MesoDyn time loop ****/
  for (int t = 0; t < timesteps; t++) {
    cout << "MESODYN: t = " << t << endl;

    New[0]->SolveMesodyn(loader_callback, solver_callback);

    //Calculate and add noise flux
    int i = 0;
    for (auto& all_fluxes : flux) {
      all_fluxes->J = solver_flux[i]->J;
      ++i;
    }

    noise_flux();
  
    //sanity_check();

    cout << "Order parameter: " << calculate_order_parameter() << endl;

    i = 0;
    for (auto all_components : component) {
      all_components->load_rho(solver_component[i]->rho);
      ++i;
    }

    if (t % timebetweensaves == 0) {
      write_output();
    }
  } // time loop

  std::cout << "Done." << std::endl;
  return true;
}

void Mesodyn::set_update_lists() {
  int c = 0;
  vector<int> temp;
  for (size_t i = 0; i < component_no; ++i) {
    for (size_t j = 0; j < (component_no - 1) - i; ++j) {
      temp.push_back(i);
      temp.push_back(c);
      update_plus.push_back(temp);
      temp.clear();
      ++c;
    }
  }

  c = 0;
  for (size_t j = 0; j < (component_no - 1); ++j) {
    for (size_t i = 1 + j; i < component_no; ++i) {
      temp.push_back(i);
      temp.push_back(c);
      update_minus.push_back(temp);
      temp.clear();
      ++c;
    }
  }
}

Real Mesodyn::calculate_order_parameter() {
  stl::device_vector<Real> difference(M);
  Real order_parameter{0};
  for (size_t i = 0; i < component_no - 1; ++i)
      for (size_t j = i + 1; j < component_no; ++j) {
          stl::transform(solver_component[i]->rho.begin(),
                         solver_component[i]->rho.end(),
                         solver_component[j]->rho.begin(),
                         difference.begin(),
                         order_param_functor());
          #ifdef PAR_MESODYN
          order_parameter = thrust::reduce(difference.begin(), difference.end(), order_parameter);
          #else
          order_parameter = std::accumulate(difference.begin(), difference.end(), order_parameter);
          #endif
      }
  order_parameter /= boundaryless_volume;

  return order_parameter;
}

int Mesodyn::noise_flux() {

  for (auto all_components : solver_component) {
    gaussian->generate(M);
    //Zero all alpha's
    stl::fill(all_components->alpha.begin(), all_components->alpha.end(), 0);
    gaussian->add_noise(all_components->alpha);
  }

  for (auto& all_fluxes : solver_flux) {
    all_fluxes->langevin_flux();
  }

  for (vector<int>& i : update_plus)
      solver_component[i[0]]->update_density(solver_flux[i[1]]->J);

  for (vector<int>& i : update_minus)
      solver_component[i[0]]->update_density(solver_flux[i[1]]->J, -1);

  int i = 0;
  for (auto& all_fluxes : flux) {
    stl::transform(all_fluxes->J.begin(), all_fluxes->J.end(), solver_flux[i]->J.begin(), all_fluxes->J.begin(), std::plus<Real>());
    ++i;
  }

  return 0;
}

int Mesodyn::sanity_check() {

  //Check local mass conservation
  Real sum {0};

  skip_bounds([this, sum](int x, int y, int z) mutable {
    sum = 0;
    for (auto all_components : solver_component) {
      sum += val(all_components->rho, x, y, z);
      //And check for negative numbers
      if ( val(all_components->rho, x, y, z) > 1 ||  val(all_components->rho, x, y, z) < 0)
        cerr << "CRITICAL ERROR: COMPONENT DENSITY > 1 || < 0" << endl;
    }
    if (sum > (1+New[0]->tolerance) || sum < (1-New[0]->tolerance)) {
      cerr << "CRITICAL ERROR: LOCAL DENSITIES DO NOT SUM TO 1, DIFFERENCE: " << sum-1 << endl;
      cerr << "Location: " << x << " " << y << " " << z << endl;
    }
  });

  //Check global mass conservation
  sum = 0;
  for (auto all_components : solver_component) {
    sum += all_components->theta();
  }

  if ( sum - boundaryless_volume > 0+New[0]->tolerance ) {
    cerr << "CRITICAL ERROR: MASS NOT CONSERVED! SUM THETA != M." << endl;
    cerr << "Difference: " << sum - ((MX-2)*(MY-2)*(MZ-2)) << endl;
  }

  return 0;
}

void Mesodyn::load_alpha(Real* alpha, size_t i) {
    if (i < component_no) {
      solver_component[i]->load_alpha(alpha);
    }
}

Real* Mesodyn::solve_explicit() {
  for (auto all_components : solver_component)
    all_components->update_boundaries();

  for (auto& all_fluxes : solver_flux)
    all_fluxes->langevin_flux();

  size_t n = 0;
  for (auto all_components : solver_component) {
    stl::copy(all_components->rho.begin(), all_components->rho.end(), rho.begin()+n*M);
    ++n;
  }

  #ifdef PAR_MESODYN
  return stl::raw_pointer_cast(&rho[0]);
  #else
  return &rho[0];
  #endif
}

Real* Mesodyn::solve_crank_nicolson() {
  int c = 0;

  for (auto all_components : solver_component) {
    all_components->load_rho(component[c]->rho);
    ++c;
  }

  for (auto all_components : solver_component)
    all_components->update_boundaries();

  for (auto& all_fluxes : solver_flux)
    all_fluxes->langevin_flux();

  for (vector<int>& i : update_plus)
      solver_component[i[0]]->update_density(flux[i[1]]->J, solver_flux[i[1]]->J, cn_ratio);


  for (vector<int>& i : update_minus)
      solver_component[i[0]]->update_density(flux[i[1]]->J, solver_flux[i[1]]->J, cn_ratio, -1);

  size_t n = 0;
  for (auto all_components : solver_component) {
    stl::copy(all_components->rho.begin(), all_components->rho.end(), rho.begin()+n*M);
    ++n;
  }

  #ifdef PAR_MESODYN
  return stl::raw_pointer_cast(&rho[0]);
  #else
  return &rho[0];
  #endif
}

int Mesodyn::initial_conditions() {

  //If molecules are pinned they cannot move, so we have to free them before moving them by using fluxes
  for (size_t i = 0; i < Seg.size(); ++i) {
    if (Seg[i]->freedom == "pinned")
      Seg[i]->freedom = "free";
  }


  stl::host_vector<stl::host_vector<Real>> rho(component_no, stl::host_vector<Real>(M));

  stl::device_vector<int> temp_mask(M);

  #if defined(PAR_MESODYN) || ! defined(CUDA)
  stl::copy(Sys[0]->KSAM, Sys[0]->KSAM+M, temp_mask.begin());
  #else
  TransferDataToHost(&temp_mask[0], Sys[0]->KSAM, M);
  #endif
  

  stl::host_vector<int> mask = temp_mask;

  switch (initialization_mode) {
    case INIT_HOMOGENEOUS:
      init_rho_homogeneous(rho, mask);
      break;
    case INIT_FROMPRO:
      {
        Reader reader(read_filename);

        if (reader.filetype != Reader::PRO)
          throw ERROR_FILE_FORMAT;

        if ( ((reader.MX)*(reader.MY)*(reader.MZ)) != (size_t)M)
          throw ERROR_SIZE_INCOMPATIBLE;

        rho = reader.multicomponent_rho;
        break;
      }
    case INIT_FROMVTK:
      for (size_t i = 0 ; i < rho.size() ; ++i) {
        string filename = read_filename + to_string(i+1) + ".vtk";
        Reader reader(filename);

        if (reader.filetype != Reader::VTK)
          throw ERROR_FILE_FORMAT;

        if ( (reader.MX)*(reader.MY)*(reader.MZ) != (size_t)M)
          throw ERROR_SIZE_INCOMPATIBLE;

        rho[i] = reader.rho;
      }
      break;
  }

  // BC0: bX0, BC1: bXm, etc.
  vector<Boundary1D::boundary> boundaries;
  for (string& boundary : Lat[0]->BC) {
    if (boundary == "mirror") {
      boundaries.push_back(Boundary1D::MIRROR);
    }
    if (boundary == "periodic") {
      boundaries.push_back(Boundary1D::PERIODIC);
    }
  }

  //Because noise per site messes up mirror boundaries, we need to tell mesodyn to mask them.
  if (boundaries[0] == Boundary1D::MIRROR) {
    x0_boundary( [this, &mask] (int x, int y, int z) mutable {
      *val_ptr(mask, x, y, z) = 0;
    });
  }
  if (boundaries[1] == Boundary1D::MIRROR) {
    xm_boundary( [this, &mask] (int x, int y, int z) mutable {
      *val_ptr(mask, x, y, z) = 0;
    });
  }
  if (dimensions > 1) {
    if (boundaries[2] == Boundary1D::MIRROR) {
      y0_boundary( [this, &mask] (int x, int y, int z) mutable {
        *val_ptr(mask, x, y, z) = 0;
      });
    }

    if (boundaries[3] == Boundary1D::MIRROR) {
      ym_boundary( [this, &mask] (int x, int y, int z) mutable {
        *val_ptr(mask, x, y, z) = 0;
      });
    }
  }
  if (dimensions > 2) {
    if (boundaries[4] == Boundary1D::MIRROR) {
      z0_boundary( [this, &mask] (int x, int y, int z) mutable {
        *val_ptr(mask, x, y, z) = 0;
      });
    }

    if (boundaries[5] == Boundary1D::MIRROR) {
      zm_boundary( [this, &mask] (int x, int y, int z) mutable {
        *val_ptr(mask, x, y, z) = 0;
      });
    }
  }

  switch (dimensions) {
  case 1:
    Mesodyn::boundary = make_shared<Boundary1D>(Lat[0], boundaries[0], boundaries[1]);

    for (size_t i = 0; i < component_no; ++i) {
      component.push_back(make_shared<Component>(Lat[0], boundary, rho[i]));
      solver_component.push_back(make_shared<Component>(Lat[0], boundary, rho[i]));
    }

    if (seed_specified == true) {
      Mesodyn::gaussian = make_unique<Gaussian_noise>(boundary, D, M, mean, stddev, seed);
    } else {
      Mesodyn::gaussian = make_unique<Gaussian_noise>(boundary, D, M, mean, stddev);
    }

    for (size_t i = 0; i < component_no - 1; ++i) {
      for (size_t j = i + 1; j < component_no; ++j) {
        flux.push_back(make_unique<Flux1D>(Lat[0], D * dt, mask, component[i], component[j]));
        solver_flux.push_back(make_unique<Flux1D>(Lat[0], D * dt, mask, solver_component[i], solver_component[j]));
      }
    }
    break;
  case 2:
    Mesodyn::boundary = make_shared<Boundary2D>(Lat[0], boundaries[0], boundaries[1], boundaries[2], boundaries[3]);

    for (size_t i = 0; i < component_no; ++i) {
      component.push_back(make_shared<Component>(Lat[0], boundary, rho[i]));
      solver_component.push_back(make_shared<Component>(Lat[0], boundary, rho[i]));
    }

    if (seed_specified == true) {
      Mesodyn::gaussian = make_unique<Gaussian_noise>(boundary, D, M, mean, stddev, seed);
    } else {
      Mesodyn::gaussian = make_unique<Gaussian_noise>(boundary, D, M, mean, stddev);
    }

    for (size_t i = 0; i < component_no - 1; ++i) {
      for (size_t j = i + 1; j < component_no; ++j) {
        flux.push_back(make_unique<Flux2D>(Lat[0], D * dt, mask, component[i], component[j]));
        solver_flux.push_back(make_unique<Flux2D>(Lat[0], D * dt, mask, solver_component[i], solver_component[j]));
      }
    }
    break;
  case 3:
      Mesodyn::boundary = make_shared<Boundary3D>(Lat[0], boundaries[0], boundaries[1], boundaries[2], boundaries[3], boundaries[4], boundaries[5]);

    for (size_t i = 0; i < component_no; ++i) {
      component.push_back(make_shared<Component>(Lat[0], boundary, rho[i]));
      solver_component.push_back(make_shared<Component>(Lat[0], boundary, rho[i]));
    }

    if (seed_specified == true) {
      Mesodyn::gaussian = make_unique<Gaussian_noise>(boundary, D, M, mean, stddev, seed);
    } else {
      Mesodyn::gaussian = make_unique<Gaussian_noise>(boundary, D, M, mean, stddev);
    }

    for (size_t i = 0; i < component_no - 1; ++i) {
      for (size_t j = i + 1; j < component_no; ++j) {
        flux.push_back(make_unique<Flux3D>(Lat[0], D * dt, mask, component[i], component[j]));
        solver_flux.push_back(make_unique<Flux3D>(Lat[0], D * dt, mask, solver_component[i], solver_component[j]));
      }
    }
    break;
  }

  norm_theta(component);
  norm_theta(solver_component);

  return 0;
}

void Mesodyn::explicit_start() {
  // Prepare callbacks
  auto explicit_solver_callback = bind(&Mesodyn::solve_explicit, this);
  auto loader_callback = bind(&Mesodyn::load_alpha, this, std::placeholders::_1, std::placeholders::_2);
  
  // start Crank-Nicolson using one explicit step
  New[0]->SolveMesodyn(loader_callback, explicit_solver_callback);

  // Update densities after first timestep
  for (vector<int>& i : update_plus)
      solver_component[i[0]]->update_density(solver_flux[i[1]]->J);

  for (vector<int>& i : update_minus)
      solver_component[i[0]]->update_density(solver_flux[i[1]]->J, -1);

  // Load rho_k+1 into rho_k
  int i = 0;
  for (auto all_components : component) {
    all_components->load_rho(solver_component[i]->rho);
    ++i;
  }

  // Load J_k+1 into J_k
  i = 0;
  for (auto& all_fluxes : flux) {
    all_fluxes->J = solver_flux[i]->J;
    ++i;
  }
}

int Mesodyn::norm_theta(vector< shared_ptr<Component> >& component) {
   size_t solvent = (size_t)Sys[0]->solvent;

  Real sum_theta{0};
  int c{0};

  vector<int> solvent_mons;

  for (size_t i = 0; i < Mol.size(); ++i) {

    size_t mon_nr = Mol[i]->MolMonList.size();

    if (i != solvent) {
      Real theta = Mol[i]->theta;
      sum_theta += theta;
      for (size_t j = 0; j < mon_nr; ++j) {

        Real mon_theta{0};
        mon_theta = theta * Mol[i]->fraction(Mol[i]->MolMonList[j]);

        Real sum_of_elements{0};
        skip_bounds([this, &sum_of_elements, component, c](int x, int y, int z) mutable {
          sum_of_elements += val(component[c]->rho, x, y, z);
        });

        skip_bounds([this, &component, c, mon_theta, sum_of_elements](int x, int y, int z) mutable {
          *val_ptr(component[c]->rho, x, y, z) *= mon_theta / sum_of_elements;
        });

        ++c;
      }
    } else {
      for (size_t j = 0; j < mon_nr; ++j) {
        solvent_mons.emplace_back(c);
        ++c;
      }
    }
  }

  // We now know the total density and can adjust the solvent accodingly to add up to 1.
  // Let's first find out how much there is to adjust.
  vector<Real> residuals(M);

  //Pool densities per position
  for (int i = 0; i < M; ++i) {
    for (size_t j = 0; j < component_no; ++j)
      residuals[i] += component[j]->rho[i];
  }

  // Calculate excesss / defecit
  for (Real& i : residuals) {
    i -= 1;
  }

  // If there's only one solvent mon, this problem is easy.
  if (solvent_mons.size() == 1) {
    skip_bounds([this, &component, residuals, solvent_mons](int x, int y, int z) mutable {
      *val_ptr(component[solvent_mons[0]]->rho, x, y, z) -= val(residuals, x, y, z);
    });
  } else {
    cerr << "Norming solvents with mutliple monomers is not supported! Please write your own script" << endl;
    throw ERROR_FILE_FORMAT;
  }

return 0;
}

int Mesodyn::init_rho_homogeneous(stl::host_vector<stl::host_vector<Real>>& rho ,stl::host_vector<int>& mask) {
  size_t solvent = (size_t)Sys[0]->solvent; // Find which mol is the solvent

  Real sum_theta{0};
  vector<Real> solvent_mons;
  int c{0};
  Real theta{0};
  Real solvent_mol{0};

  for (size_t i = 0; i < Mol.size(); ++i) {
    size_t mon_nr = Mol[i]->MolMonList.size();
    if (i != solvent) {
      theta = Mol[i]->theta;
      sum_theta += theta;
      for (size_t j = 0; j < mon_nr; ++j) {
        Real mon_theta{0};
        mon_theta = theta * Mol[i]->fraction(Mol[i]->MolMonList[j]);
        for (int z = 0; z < M; ++z) {
          for (size_t i = 0; i < component_no; ++i) {
            if (mask[z] == 0)
              rho[i][z] = 0;
            else
              rho[c][z] = mon_theta / boundaryless_volume;
          }
        }
        ++c;
      }
    } else {
      solvent_mol = i;
      for (size_t j = 0; j < mon_nr; ++j) {
        solvent_mons.push_back(c);
        ++c;
      }
    }
  }

  for (size_t i = 0; i < solvent_mons.size(); ++i) {
    Real mon_theta{0};
    mon_theta = (boundaryless_volume - sum_theta) * Mol[solvent_mol]->fraction(Mol[solvent_mol]->MolMonList[i]);
    for (int z = 0; z < M; ++z) {
      if (mask[z] == 0)
        rho[solvent][z] = 0;
      else
        rho[solvent_mons[i]][z] = mon_theta / boundaryless_volume;
    }
  }

  return 0;
}

/******* Output generation *******/

void Mesodyn::set_filename() {
  /* Open filestream and set filename to "mesodyn-datetime.csv" */
  string output_folder = "output/";
  string bin_folder = "bin";

  // Find path to Namics executable
  char result[PATH_MAX];
  ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
  string executable_path = string(result, (count > 0) ? count : 0);

  // Find the last string before the executable
  size_t found = executable_path.find_last_of("/\\");

  // Set the output folder to be one level up from the binary folder, plus the specified output folder
  output_folder = executable_path.substr(0, found - bin_folder.size()) + output_folder;

  filename << output_folder << "mesodyn-";

  time_t rawtime;
  time(&rawtime);
  filename << rawtime;
}

int Mesodyn::write_output() {
  //Implement below for more starts? (OutputList higher numbers?)

    if (!Out[0]->CheckInput(1)) {
      cout << "input_error in output " << endl;
      return 0;
    }

    Out[0]->output_nr = writes;
    Out[0]->n_output = timesteps / timebetweensaves;

    New[0]->PushOutput();

    Out[0]->push("order_parameter",order_parameter);
    Out[0]->push("timesteps", timesteps);
    Out[0]->push("timebetweensaves", timebetweensaves);
    Out[0]->push("diffusionconstant", D);
    Out[0]->push("seed", seed);
    Out[0]->push("mean", mean);
    Out[0]->push("stddev", stddev);
    Out[0]->push("delta_t", dt);
    Out[0]->push("cn_ratio", cn_ratio);

    int component_count{1};
    for (auto all_components : solver_component) {
      stl::host_vector<Real> t_rho = all_components->rho;
      ostringstream vtk_filename;
      vtk_filename << filename.str() << "-rho" << component_count << "-" << writes << ".vtk";
      Out[0]->vtk_structured_grid(vtk_filename.str(), &t_rho[0], component_count);
      ++component_count;
    }
    ++writes;

    Out[0]->WriteOutput(writes);

  return 0;
}

/******* FLUX: TOOLS FOR CALCULATING FLUXES BETWEEN 1 PAIR OF COMPONENTS, HANDLING OF SOLIDS *********/

Flux1D::Flux1D(Lattice* Lat, Real D, stl::host_vector<int>& mask, shared_ptr<Component> A, shared_ptr<Component> B)
    : Lattice_Access(Lat), J_plus(M), J_minus(M), J(M), A{A}, B{B}, L(M), mu(M), D{D}, JX{Lat->JX} {
  Flux1D::mask(mask);
}

Flux2D::Flux2D(Lattice* Lat, Real D, stl::host_vector<int>& mask, shared_ptr<Component> A, shared_ptr<Component> B)
    : Flux1D(Lat, D, mask, A, B), JY{Lat->JY} {
  Flux2D::mask(mask);
}

Flux3D::Flux3D(Lattice* Lat, Real D, stl::host_vector<int>& mask, shared_ptr<Component> A, shared_ptr<Component> B)
    : Flux2D(Lat, D, mask, A, B), JZ{Lat->JZ} {
  Flux3D::mask(mask);
}

Flux1D::~Flux1D() {
}

Flux2D::~Flux2D() {
}

Flux3D::~Flux3D() {
}

Real Flux1D::J_at(int x, int y, int z) {
  return val(J, x, y, z);
}

Real Flux1D::L_at(int x, int y, int z) {
  return val(L, x, y, z);
}

Real Flux1D::mu_at(int x, int y, int z) {
  return val(mu, x, y, z);
}

int Flux1D::mask(stl::host_vector<int>& mask_in) {
  if ((int)mask_in.size() != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  // No fluxes will ever be calculated going from the boundary into the system
  skip_bounds([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x, y, z) == 1) {
      if (val(mask_in, x + 1, y, z) == 1) {
        Mask_plus_x.push_back(index(x, y, z));
      }
      if (val(mask_in, x - 1, y, z) == 1) {
        Mask_minus_x.push_back(index(x, y, z));
      }
    }
  });

  //for the x-boundary:
  x0_boundary([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x+1, y, z) == 1 && val(mask_in, x, y, z) == 1)
    Mask_plus_x.push_back(index(x, y, z));
  });
  xm_boundary([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x-1, y, z) == 1 && val(mask_in, x, y, z) == 1)
    Mask_minus_x.push_back(index(x, y, z));
  });

  return 0;
}

int Flux2D::mask(stl::host_vector<int>& mask_in) {
  if ((int)mask_in.size() != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  // No fluxes will ever be calculated going from the boundary into the system
  skip_bounds([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x, y, z) == 1) {
      if (val(mask_in, x, y + 1, z) == 1)
        Mask_plus_y.push_back(index(x, y, z));
      if (val(mask_in, x, y - 1, z) == 1)
        Mask_minus_y.push_back(index(x, y, z));
    }
  });

  //for the y-boundary:
  y0_boundary([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x, y+1, z) == 1 && val(mask_in, x, y, z) == 1)
    Mask_plus_y.push_back(index(x, y, z));
  });
  ym_boundary([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x, y-1, z) == 1 && val(mask_in, x, y, z) == 1)
    Mask_minus_y.push_back(index(x, y, z));
  });

  return 0;
}

int Flux3D::mask(stl::host_vector<int>& mask_in) {
  if ((int)mask_in.size() != M) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  // No fluxes will ever be calculated going from the boundary into the system, messing up periodic boundaries
  skip_bounds([this, mask_in](int x, int y, int z) mutable {
    if (val(mask_in, x, y, z) == 1) {
      if (val(mask_in, x, y, z + 1) == 1) {
        Mask_plus_z.push_back(index(x, y, z));
      }
      if (val(mask_in, x, y, z - 1) == 1) {
        Mask_minus_z.push_back(index(x, y, z));
      }
    }
  });

  //for the z-boundary:
  z0_boundary([this, mask_in](int x, int y, int z) mutable {
      if (val(mask_in, x, y, z+1) == 1 && val(mask_in, x, y, z) == 1)
        Mask_plus_z.push_back(index(x, y, z));
  });
  zm_boundary([this, mask_in](int x, int y, int z) mutable {
      if (val(mask_in, x, y, z-1) == 1 && val(mask_in, x, y, z) == 1)
        Mask_minus_z.push_back(index(x, y, z));
  });

  return 0;
}

int Flux1D::langevin_flux() {

  //Zero (with bounds checking) vector J before use
  stl::fill(J.begin(), J.end(), 0);

  if (A->rho.size() != J.size()) {
    //already checked: A.alpha.size = B.alpha.size and A.rho.size = B.rho.size
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  onsager_coefficient(A->rho, B->rho);
  potential_difference(A->alpha, B->alpha);
  langevin_flux(Mask_plus_x, Mask_minus_x, JX);

  return 0;
}

int Flux2D::langevin_flux() {
  Flux1D::langevin_flux();

  Flux1D::langevin_flux(Mask_plus_y, Mask_minus_y, JY);

  return 0;
}

int Flux3D::langevin_flux() {
  Flux2D::langevin_flux();

  Flux1D::langevin_flux(Mask_plus_z, Mask_minus_z, JZ);

  return 0;
}

int Flux1D::onsager_coefficient(stl::device_vector<Real>& A, stl::device_vector<Real>& B) {

  if (A.size() != B.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }
  stl::transform(A.begin(), A.end(), B.begin(), L.begin(), stl::multiplies<Real>());

  return 0;
}

int Flux1D::potential_difference(stl::device_vector<Real>& A, stl::device_vector<Real>& B) {

  if (A.size() != B.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  stl::transform(A.begin(), A.end(), B.begin(), mu.begin(), stl::minus<Real>());

  return 0;
}

int Flux1D::langevin_flux(stl::host_vector<int>& mask_plus, stl::host_vector<int>& mask_minus, int jump) {

  //////
  //
  // All these transform and fill functions are needed to make mesodyn parallelizable using Thrust
  //
  // The code below does the equivalent of the following code (and faster):
  //  
  //   for (vector<int>::iterator itt = mask_plus.begin() ; itt < mask_plus.end(); ++itt) {
  //     auto z = *itt;
  //     J_plus[z] = -D * ((L[z] + L[z + jump]) * (mu[z + jump] - mu[z]));
  //  }
  //
  //  for (vector<int>::iterator itt = mask_minus.begin() ; itt < mask_minus.end(); ++itt) {
  //     auto z = *itt;
  //     J_minus[z] = -J_plus[z - jump];
  //  }
  //
  //////

  stl::fill(J_minus.begin(), J_minus.end(), 0);
  stl::fill(J_plus.begin(), J_plus.end(), 0);

  stl::device_vector<Real> t_mu(M);
  stl::device_vector<Real> t_L(M);

  stl::transform(mu.begin()+jump, mu.end(), mu.begin(), t_mu.begin(), stl::minus<Real>());
  stl::transform(L.begin(), L.end()-jump, L.begin()+jump, t_L.begin(), stl::plus<Real>());
  stl::transform(t_mu.begin(), t_mu.end(), t_L.begin(), J_plus.begin(), const_multiply_functor(-D) );

  stl::transform(J_plus.begin(), J_plus.end()-jump, J_minus.begin()+jump, stl::negate<Real>());

  stl::transform(J_plus.begin(), J_plus.end(), J.begin(), J.begin(), stl::plus<Real>());
  stl::transform(J_minus.begin(), J_minus.end(), J.begin(), J.begin(), stl::plus<Real>());

  return 0;
}

/****************** COMPONENT: DENSITY PROFILE STORAGE AND UPDATING, BOUNDARY CONDITIONS ********************/

/******* Constructors *******/

Component::Component(Lattice* Lat, shared_ptr<Boundary1D> boundary, stl::host_vector<Real>& rho)
    : Lattice_Access(Lat), rho(M), alpha(M), boundary(boundary) {
  //This check is implemented multiple times throughout mesodyn because rho and alpha are public.
  if (rho.size() != alpha.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  this->rho = rho;

  #ifdef PAR_MESODYN
  rho_ptr = thrust::raw_pointer_cast(&this->rho[0]);
  alpha_ptr = thrust::raw_pointer_cast(&this->alpha[0]);
  #else
  rho_ptr = &Component::rho[0];
  alpha_ptr = &Component::alpha[0];
  #endif

  update_boundaries();
}

Component::~Component() {
}

/******* Interface *******/

Real Component::rho_at(int x, int y, int z) {
  return val(rho, x, y, z);
}

Real Component::alpha_at(int x, int y, int z) {
  return val(alpha, x, y, z);
}

int Component::update_density(stl::device_vector<Real>& J, int sign) {
  //Explicit update

  if (J.size() != rho.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }

  stl::transform(J.begin(), J.end(), rho.begin(), rho.begin(), saxpy_functor(sign) );

  return 0;
}

int Component::update_density(stl::device_vector<Real>& J1, stl::device_vector<Real>& J2, Real ratio, int sign) {
  //Implicit update
  if (J1.size() != rho.size() || J1.size() != J2.size()) {
    throw ERROR_SIZE_INCOMPATIBLE;
  }
  // Rho <- A * J1 + Rho
  stl::transform(J1.begin(), J1.end(), rho.begin(), rho.begin(), saxpy_functor(sign*ratio) );
  stl::transform(J2.begin(), J2.end(), rho.begin(), rho.begin(), saxpy_functor((1.0f-ratio)*sign) );

  return 0;
}

int Component::load_alpha(stl::device_vector<Real> alpha) {
  stl::copy(alpha.begin(), alpha.end(), Component::alpha.begin());
  return 0;
}

int Component::load_alpha(Real* alpha) {
  #if defined(CUDA) && ! defined(PAR_MESODYN)
  TransferDataToHost(&Component::alpha[0], alpha, M);
  #else
  stl::copy(alpha, alpha+M, this->alpha.begin());
  #endif
  return 0;
}

int Component::load_rho(stl::device_vector<Real> rho) {
  stl::copy(rho.begin(), rho.end(), this->rho.begin());
  return 0;
}

int Component::load_rho(Real* rho) {
  #if defined(CUDA) && ! defined(PAR_MESODYN)
  TransferDataToHost(&Component::rho[0], rho, M);
  #else
  stl::copy(rho, rho+M, this->rho.begin());
  #endif
  return 0;
}

int Component::update_boundaries() {
  #ifdef PAR_MESODYN
  Lat->set_bounds(alpha_ptr);
  Lat->set_bounds(rho_ptr);
  #else
  boundary->update_boundaries(alpha);
  boundary->update_boundaries(rho);
  #endif

  return 0;
}

Real Component::theta() {
  Real sum{0};
  #ifdef PAR_MESODYN
  Lat->remove_bounds(rho_ptr);
  sum = stl::reduce(rho.begin(), rho.end(),0);
  #else
  skip_bounds([this, &sum](int x, int y, int z) mutable {
    sum += val(rho, x, y, z);
  });
  #endif
  return sum;
}

/****************** Lattice_Access: AN INTERFACE FOR LATTICE ********************/

Lattice_Access::Lattice_Access(Lattice* Lat)
    : dimensions{Lat->gradients}, Lat{Lat}, JX{Lat->JX}, JY{Lat->JY}, JZ{Lat->JZ}, M{Lat->M}, MX{2 + Lat->MX}, MY{setMY(Lat)}, MZ{setMZ(Lat)} {
  // If this is not true, NOTHING will work. So this check is aboslutely necessary.
  // If, at some point, the above happens to be the case every class in this module will probably need rewriting.
  assert(
      (MX > 0 && MY == 0 && MZ == 0) ||
      (MX > 0 && MY > 0 && MZ == 0) ||
      (MX > 0 && MY > 0 && MZ > 0));
}

Lattice_Access::~Lattice_Access() {
}

inline void Lattice_Access::par_skip_bounds(function<void(int, int, int)> function) {
  int x{1};
  int y{1};

  #pragma omp parallel for private(y,x)
  for ( int z = 1 ; z < MZ-1 ; ++z ) {
      y = 1;
      do {
        x = 1;
        do {
          function(x, y, z);
          ++x;
        } while (x < MX - 1);
        ++y;
      } while (y < MY - 1);
  }
}

inline void Lattice_Access::skip_bounds(function<void(int, int, int)> function) {
  int x{1};
  int y{1};
  int z{1};

  do {
      y = 1;
      do {
        x = 1;
        do {
          function(x, y, z);
          ++x;
        } while (x < MX - 1);
        ++y;
      } while (y < MY - 1);
      ++z;
    } while (z < MZ - 1);
}

inline void Lattice_Access::bounds(function<void(int, int, int)> function) {

  x0_boundary(function);
  xm_boundary(function);

  if (dimensions > 1) {
    y0_boundary(function);
    ym_boundary(function);
  }

  if (dimensions > 2) {
    z0_boundary(function);
    zm_boundary(function);
  }
}

inline void Lattice_Access::x0_boundary(function<void(int, int, int)> function) {
  int x = 0;
  int y = 0;
  int z = 0;
  do {
    y = 0;
    do {
      function(x, y, z);
      ++y;
    } while (y < MY);
    ++z;
  } while (z < MZ);
}

inline void Lattice_Access::xm_boundary(function<void(int, int, int)> function) {
  int x = MX - 1, y = 0, z = 0;
  do {
    y = 0;
    do {
      function(x, y, z);
      ++y;
    } while (y < MY);
    ++z;
  } while (z < MZ);
}

inline void Lattice_Access::y0_boundary(function<void(int, int, int)> function) {
  int x = 0, y = 0, z = 0;
  do {
    x = 0;
    do {
      function(x, y, z);
      ++x;
    } while (x < MX);
    ++z;
  } while (z < MZ);
}

inline void Lattice_Access::ym_boundary(function<void(int, int, int)> function) {
  int x = 0, y = MY - 1, z = 0;
  do {
    x = 0;
    do {
      function(x, y, z);
      ++x;
    } while (x < MX);
    ++z;
  } while (z < MZ);
}

inline void Lattice_Access::z0_boundary(function<void(int, int, int)> function) {
  int x = 0, y = 0, z = 0;
  do {
    x = 0;
    do {
      function(x, y, z);
      ++x;
    } while (x < MX);
    ++y;
  } while (y < MY);
}

inline void Lattice_Access::zm_boundary(function<void(int, int, int)> function) {
  int x = 0, y = 0, z = MZ - 1;
  do {
    x = 0;
    do {
      function(x, y, z);
      ++x;
    } while (x < MX);
    ++y;
  } while (y < MY);
}

int Lattice_Access::setMY(Lattice* Lat) {
  //used by constructor
  if (dimensions < 2)
    return 0;
  else
    return Lat->MY + 2;
}

int Lattice_Access::setMZ(Lattice* Lat) {
  //used by constructor
  if (dimensions < 3)
    return 0;
  else
    return Lat->MZ + 2;
}

inline int Lattice_Access::index(int x, int y, int z) {
  return x * JX + y * JY + z * JZ;
}

inline vector<int> Lattice_Access::coordinate(int n) {
  int x = 0, y = 0, z = 0;
  int mod = 0;
  x = n / JX;
  mod = n % JX;

  if (dimensions > 1) {
    y = mod / JY;
    mod = mod % JY;
  }

  if (dimensions > 2) {
    z = mod / JZ;
  }
  return {x, y, z};
}

/******* Boundary conditions *******/

Boundary1D::Boundary1D(Lattice* Lat, boundary x0, boundary xm)
    : Lattice_Access(Lat) {
  set_x_boundaries(x0, xm);
}

Boundary2D::Boundary2D(Lattice* Lat, boundary x0, boundary xm, boundary y0, boundary ym)
    : Boundary1D(Lat, x0, xm) {
  set_y_boundaries(y0, ym);
}

Boundary3D::Boundary3D(Lattice* Lat, boundary x0, boundary xm, boundary y0, boundary ym, boundary z0, boundary zm)
    : Boundary2D(Lat, x0, xm, y0, ym) {
  set_z_boundaries(z0, zm);
}

Boundary1D::~Boundary1D() {}

Boundary2D::~Boundary2D() {}

Boundary3D::~Boundary3D() {}

int Boundary1D::update_boundaries(vector<Real>& target) {
  bX0(target);
  bXm(target);
  return 0;
}

int Boundary2D::update_boundaries(vector<Real>& target) {
  Boundary1D::update_boundaries(target);
  bY0(target);
  bYm(target);
  return 0;
}

int Boundary3D::update_boundaries(vector<Real>& target) {
  Boundary2D::update_boundaries(target);
  bZ0(target);
  bZm(target);

  return 0;
}

int Boundary1D::set_x_boundaries(boundary x0, boundary xm, Real bulk) {
  using namespace std::placeholders;

  switch (x0) {
  case MIRROR:
    bX0 = bind(&Boundary1D::bX0Mirror, this, _1);
    break;
  case PERIODIC:
    if (x0 != xm) {
      throw ERROR_PERIODIC_BOUNDARY;
    }
    bX0 = bind(&Boundary1D::bXPeriodic, this, _1);
    break;
  }

  switch (xm) {
  case MIRROR:
    bXm = bind(&Boundary1D::bXmMirror, this, _1);
    break;
  case PERIODIC:
    if (x0 != xm) {
      throw ERROR_PERIODIC_BOUNDARY;
    }
    //TODO: fix this double periodic (including ones below)
    bXm = bind(&Boundary1D::bXPeriodic, this, _1);
    break;
  }

  return 0;
}

int Boundary2D::set_y_boundaries(boundary y0, boundary ym, Real bulk) {
  using namespace std::placeholders;

  switch (y0) {
  case MIRROR:
    bY0 = bind(&Boundary2D::bY0Mirror, this, _1);
    break;
  case PERIODIC:
    if (y0 != ym) {
      throw ERROR_PERIODIC_BOUNDARY;
    }
    bY0 = bind(&Boundary2D::bYPeriodic, this, _1);
    break;
  }

  switch (ym) {
  case MIRROR:
    bYm = bind(&Boundary2D::bYmMirror, this, _1);
    break;
  case PERIODIC:
    if (y0 != ym) {
      throw ERROR_PERIODIC_BOUNDARY;
    }
    bYm = bind(&Boundary2D::bYPeriodic, this, _1);
    break;
  }

  return 0;
}

int Boundary3D::set_z_boundaries(boundary z0, boundary zm, Real bulk) {
  using namespace std::placeholders;

  switch (z0) {
  case MIRROR:
    bZ0 = bind(&Boundary3D::bZ0Mirror, this, _1);
    break;
  case PERIODIC:
    if (z0 != zm) {
      throw ERROR_PERIODIC_BOUNDARY;
    }
    bZ0 = bind(&Boundary3D::bZPeriodic, this, _1);
    break;
  }

  switch (zm) {
  case MIRROR:
    bZm = bind(&Boundary3D::bZmMirror, this, _1);
    break;
  case PERIODIC:
    if (z0 != zm) {
      throw ERROR_PERIODIC_BOUNDARY;
    }
    bZm = bind(&Boundary3D::bZPeriodic, this, _1);
    break;
  }

  return 0;
}

void Boundary1D::bX0Mirror(vector<Real>& target) {
  x0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x + 1, y, z); //start
  });
}

void Boundary1D::bXmMirror(vector<Real>& target) {
  xm_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x - 1, y, z); //end
  });
}

void Boundary1D::bXPeriodic(vector<Real>& target) {

  x0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, MX - 2, y, z); //start
  });

  xm_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, 1, y, z); //end
  });

}

void Boundary2D::bY0Mirror(vector<Real>& target) {
  y0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y + 1, z); //start
  });
}

void Boundary2D::bYmMirror(vector<Real>& target) {
  ym_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y - 1, z); //end
  });
}

void Boundary2D::bYPeriodic(vector<Real>& target) {
  y0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, MY - 2, z); //start
  });
  ym_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, 1, z); //end
  });
}
void Boundary3D::bZ0Mirror(vector<Real>& target) {
  z0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y, z + 1); //start
  });
}

void Boundary3D::bZmMirror(vector<Real>& target) {
  zm_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y, z - 1); //end
  });
}

void Boundary3D::bZPeriodic(vector<Real>& target) {
  z0_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y, MZ - 2); //start
  });
  zm_boundary([this, &target](int x, int y, int z) mutable {
    *val_ptr(target, x, y, z) = val(target, x, y, 1); //end
  });
}

/******* GAUSSIAN_NOISE: GENERATES WHITE NOISE FOR FLUXES ********/

Gaussian_noise::Gaussian_noise(shared_ptr<Boundary1D> boundary, Real D, int size, Real mean, Real stddev) : noise(size), prng{std::random_device{}()}, dist(mean, stddev), boundary{boundary} {}

Gaussian_noise::Gaussian_noise(shared_ptr<Boundary1D> boundary, Real D, int size, Real mean, Real stddev, size_t seed) : noise(size), prng(seed), dist(mean, stddev), boundary{boundary} {}

int Gaussian_noise::generate(size_t M) {
  stl::host_vector<Real> tmp_noise(M);
  for (Real& value : tmp_noise)
    value = dist(prng);

  noise = tmp_noise;

  #ifdef PAR_MESODYN
  boundary->Lat->set_bounds( thrust::raw_pointer_cast(&noise[0]) );
  #else
  boundary->update_boundaries(noise);
  #endif

  return 0;
}

int Gaussian_noise::add_noise(stl::device_vector<Real>& target) {
  stl::transform(noise.begin(), noise.end(), target.begin(), target.begin(), stl::plus<Real>());
  return 0;
}

/******* TOOLS: MATHEMATICS AND INTERFACE FOR THE IO CLASSES ********/

/***** Mathematics *****/

int Mesodyn::factorial(int n) {
  if (n > 1) {
    return n * factorial(n - 1);
  } else
    return 1;
}

int Mesodyn::combinations(int n, int k) {
  return factorial(n) / (factorial(n - k) * factorial(k));
}

/***** IO *****/

int Mesodyn::GetValue(string prop, int& int_result, Real& Real_result, string& string_result) {
  int i = 0;
  int length = ints.size();
  while (i < length) {
    if (prop == ints[i]) {
      int_result = ints_value[i];
      return 1;
    }
    i++;
  }
  i = 0;
  length = Reals.size();
  while (i < length) {
    if (prop == Reals[i]) {
      Real_result = Reals_value[i];
      return 2;
    }
    i++;
  }
  i = 0;
  length = bools.size();
  while (i < length) {
    if (prop == bools[i]) {
      if (bools_value[i])
        string_result = "true";
      else
        string_result = "false";
      return 3;
    }
    i++;
  }
  i = 0;
  length = strings.size();
  while (i < length) {
    if (prop == strings[i]) {
      string_result = strings_value[i];
      return 3;
    }
    i++;
  }
  return 0;
}

string Mesodyn::GetValue(string parameter){
	int i=0;
	int length = PARAMETERS.size();
	while (i<length) {
		if (parameter==PARAMETERS[i]) {
			return VALUES[i];
		}
		i++;
	}
	return "" ;
}
