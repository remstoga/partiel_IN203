# Examen Machine

Rémi GASTAUD

Nombre de coeurs physiques : 6

Nombre de threads : 12

Toutes les mesures sont effectuées avec les valeurs par défaut

## Mesure du temps en séquentiel

| DFT          | Sélection des p% coefficients | Reconstitution de l'image |
| ------------ | ----------------------------- | ------------------------- |
| 102 secondes | 5 ms                          | 12 secondes               |

C'est la première DFT qui prend énormément de temps, ensuite vient la reconstitution de l'image, et enfin la sélection des coefficients est négligeable. C'est un résultat peu surprenants puisque la reconstitution revient quasiment au même en temps de calcul que la DFT, sauf qu'elle n'est effectuée que sur 10% des coefficients, donc 10 fois mois longue.

## Parallélisation OpenMP

On choisit donc de paralléliser les boucles de la DFT et de la reconstitution.

Pour la boucle DFT, on parallélise sur les pixels, et pour celle de la reconstitution, sur les fréquences. Cela permet de ne pas avoir de problèmes avec des conflits d'accès mémoire.

On obtient les résultats suivants :

| Nombre de threads | DFT         | Reconstitution de l'image |
| ----------------- | ----------- | ------------------------- |
| 2                 | 61 secondes | 6.5 secondes              |
| 4                 | 32 secondes | 3.7 secondes              |
| 5                 | 27 secondes | 3 secondes                |
| 6                 | 23 secondes | 2.5 secondes              |
| 7                 | 25 secondes | 2.9 secondes              |
| 8                 | 25 secondes | 2.7 secondes              |
| 10                | 19 secondes | 1.9                       |
| 12                | 19 secondes | 2.3                       |

On remarque qu'à partir de 6/7 threads, les résultats ne s'améliorent que très peu : cela est du à la loi d'Amdhal, qui indique la limite de l'accélération d'un code parallélisé, et à la limite en accès mémoire. Ici on peut en déduire qu'environ les trois quarts de la DFT sont parallélisables, idem pour la reconstitution de l'image.

## Première parallélisation MPI

La première parallélisation MPI n'a d'effet que sur la première transformée de Fourier, donc nous ne mesurerons que cela.

On obtient les résultats suivants :

| Nombre de processus | DFT         |
| ------------------- | ----------- |
| 2                   | 52 secondes |
| 4                   | 28 secondes |
| 5                   | 23 secondes |
| 6                   | 21 secondes |
| 7                   | 25 secondes |
| 8                   | 22 secondes |
| 10                  | 16 secondes |
| 12                  | 15 secondes |

On observe le même phénomène qu'avec openMP, et pour les même raisons, mais MPI semble être un peu plus efficace pour gagner du temps.

à noter que si le nombre de processus ne divise pas la taille de l'image, l'image compressée aura des bandes erronées en bas.

## Deuxième parallélisation MPI

Partitionner sur les fréquences permet à chaque processus d'effectuer tous les calculs de son côté, à l'exception de la recherche des coeff maximaux. Pour cela deux choix s'offrent à nous, on peut soit partager à tous les processus les valeurs de tous les coefficients dans une étape intermédiaires, puis chacun calcule les plus grands, soit partager seulement les plus grands coeff de chaque processus, les partager, et recalculer les plus grands coeff. Je n'ai pas eu le temps de mettre en oeuvre la deuxième option, donc je reste sur la première.

On obtient les résultats suivants:

| Nombre de processus | DFT           | Reconstitution de l'image |
| ------------------- | ------------- | ------------------------- |
| 2                   | 53 secondes   | 6.8 secondes              |
| 4                   | 26 secondes   | 6.7 secondes              |
| 5                   | 23 secondes   | 6.9                       |
| 6                   | 20 secondes   | 6.3   secondes            |
| 7                   | 26 secondes   | 6.2                       |
| 8                   | 22 secondes   | 9 secondes                |
| 10                  | 18 secondes   | 6.0 secondes              |
| 12                  | 15   secondes | 6.0 secondes              |

Sans grande surprise, les résultats pour la DFT sont les mêmes qu'à la question précédente, toutefois la reconstitution de l'image n'arrive pas à aller au-delà de 6 secondes, alors qu'openMP pouvait aller bien plus loin. Cela vient possiblement du fait que l'on parallélise sur les fréquences, et non les pixels comme on avait pu le faire avec openMP : le temps demandé pour les calculs est très hétérogène entre les processus, et la réduction bloquatne MPI_Reduce nécessite d'attendre tous les processus. Ainsi c'est le plus long temps qui sera mesuré.

On remarque également un pic pour 8 processus : cela doit être du à la façon dont est agencée la mémoire cache : avec 8 processus, on accède tout le temps au même emplacement mémoire, et celu-ci doit être réagencé à chaque accès.

