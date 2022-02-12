# Plot of distribution of lifetime reproductive success.

plot_lrs_distr_allpop <- function(population.name, distr, theoretic_distri){
  if(population.name == "Hutt"){
    plot_lrs_distr_hutt(distr, theoretic_distri)
  }else if(population.name == "Ache"){
    plot_lrs_distr_ache(distr, theoretic_distri)
  }else if(population.name == "Hadza"){
    plot_lrs_distr_hadza(distr, theoretic_distri)
  }else if(population.name == "Roe deer"){
    plot_lrs_distr_roedeer(distr)
  }else if(population.name == "Tsuga"){
    plot_lrs_distr_tsuga(distr)
  }
}
