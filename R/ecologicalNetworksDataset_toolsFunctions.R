dataset_names_MA <- c('Sinohara_1_ALL_PH', 'Sinohara_2_ALL_PH', 'Sinohara_3_ALL_PH', 'Sinohara_4_ALL_PH', 'Sinohara_ALL_A_PH',
                      'Sinohara_ALL_E_PH', 'Sinohara_ALL_I_PH', 'Sinohara_2_E_PH', 'Sinohara_3_E_PH', 'Sinohara_4_I_PH',
                      'Melian_OO_OO_PH', 'Hackett_1_ALL_PH', 'Hackett_2_ALL_PH', 'Hackett_1_S_PH', 'Hackett_1_GL_PH',
                      'Pocock_OO_OO_PH','Melian_OO_OO_HSD')
dataset_names_MM <-c('Melian_OO_OO_PSD','Dattilo_OO_OO_PSD','Dattilo_OO_OO_PA')
dataset_names_AA <-c('McFayden_ALL_A_HPa', 'McFayden_1_A_HPa', 'McFayden_2_A_HPa', 'McFayden_3_A_HPa', 'McFayden_4_A_HPa',
                     'McFayden_5_A_HPa', 'McFayden_6_A_HPa', 'McFayden_7_A_HPa', 'McFayden_8_A_HPa', 'McFayden_9_A_HPa',
                     'McFayden_10_A_HPa', 'McFayden_ALL_B_HPa', 'McFayden_1_B_HPa', 'McFayden_2_B_HPa', 'McFayden_3_B_HPa',
                     'McFayden_4_B_HPa', 'McFayden_5_B_HPa', 'McFayden_6_B_HPa', 'McFayden_7_B_HPa', 'McFayden_8_B_HPa',
                     'McFayden_9_B_HPa', 'McFayden_10_B_HPa', 'Hackett_1_ALL_HPa', 'Hackett_1_WL_HPa')


#' tool 1
#'
#' @param x à faire
#' @param y à faire
#' @return à faire
#' @examples
#' à faire
#'
#' @export
find_partites_type=function(dataset_name){
  if (startsWith( dataset_name, 'Ibanez')){
    interactions = list("pollination", "herbivory")
    linking_set='plants'
  } else if (startsWith( dataset_name, 'Sinohara')){
    interactions = list("pollination", "herbivory")
    linking_set='plants'
  } else if (startsWith( dataset_name, 'Melian')){
    if (str_detect( dataset_name, 'HSD')) {
      interactions = list("dispersion", "herbivory")
      linking_set='plants'
    } else if (str_detect( dataset_name, 'PH')) {
      interactions = list("pollination", "herbivory")
      linking_set='plants'
    } else if (str_detect(dataset_name, 'PSD')) {
      interactions = list("pollination","dispersion")
      linking_set='plants'
    }
  } else if (startsWith( dataset_name, 'Hackett')){
    if (str_detect(dataset_name,'_PH')){
      interactions = list("pollination", "herbivory")
      linking_set='plants'
    } else if (str_detect(dataset_name,'_HPa')) {
      interactions = list('herbivory', 'parasitism')
      linking_set='hosts'
    } else if (str_detect(dataset_name,'_SHPa')) {
      interactions = list('herbivory', 'parasitism')
      linking_set='hosts'
    } else if (str_detect(dataset_name,'_LHSH')) {
      interactions = list('leaf_herbivory', 'seed_herbivory')
      linking_set='plants'
    }
  } else if (startsWith( dataset_name, 'Pocock')){
    if (str_detect(dataset_name,'_PH')){
      interactions = list("pollination", "herbivory")
      linking_set='plants'
    } else if (str_detect(dataset_name,'_SHH')){
      interactions = list("seed_herbivory", "herbivory")
      linking_set='plants'
    }
  } else if (startsWith( dataset_name, 'Genrich')){
    interactions = list("dispersion","seed_herbivory")
    linking_set='plants'
  } else if (startsWith( dataset_name, 'Melian_DH')){
    interactions = list("dispersion", "herbivory")
    linking_set='plants'
  } else if (startsWith( dataset_name,  'Mello')){
    interactions = list("frugivory", "nectarivory")
    linking_set='plants'
  } else if (str_detect( dataset_name,  '_PSD')) {
    interactions = list("dispersion", "pollination")
    linking_set='plants'
  } else if (str_detect( dataset_name,  '_PA')){
    interactions = list("ant", "pollination")
    linking_set='plants'
  } else if (str_detect( dataset_name,  '_SDA')){
    interactions = list("ant", "dispersion")
    linking_set='plants'
  } else if (startsWith( dataset_name,  'McFayden')){
    interactions = list('herbivory', 'parasitism')
    linking_set='hosts'
  } else if (startsWith( dataset_name,  'Dattilo')){
    if (endsWith( dataset_name, 'PSD')) {
      interactions = list("pollination","dispersion")
      linking_set='plants'
    }
    else if (str_detect( dataset_name, 'PA')){
      interactions = list("pollination","ant")
      linking_set='plants'
    }
  } else {
    interactions = list("pollination", "herbivory")
    linking_set='plants'
  }
  return(list("interactions"=interactions,"linking_set"=linking_set))
}
