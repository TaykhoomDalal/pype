#!/bin/bash

sbatch visualize.sh -abdomen not_annotated --no_annotations --plot_manhattan
sbatch visualize.sh -heart not_annotated --no_annotations --plot_manhattan
sbatch visualize.sh -eyes not_annotated --no_annotations --plot_manhattan
sbatch visualize.sh -brain not_annotated --no_annotations --plot_manhattan

sbatch visualize.sh -abdomen annotated --plot_manhattan
sbatch visualize.sh -heart annotated --plot_manhattan
sbatch visualize.sh -eyes annotated --plot_manhattan
sbatch visualize.sh -brain annotated --plot_manhattan

sbatch visualize.sh -abdomen not_annotated --no_annotations --plot_bar
sbatch visualize.sh -heart not_annotated --no_annotations --plot_bar
sbatch visualize.sh -eyes not_annotated --no_annotations --plot_bar
sbatch visualize.sh -brain not_annotated --no_annotations --plot_bar