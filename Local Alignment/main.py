import argparse
import json

import local_alignment

parser = argparse.ArgumentParser(description="Local alignment configuration")
parser.add_argument("mode", type=str, help="Execution mode: dynprog / dynproglin / heuralign")
parser.add_argument("alphabet", type=str, help="Sequence alphabet")
parser.add_argument("matrix", type=str, help="Substitution matrix")
parser.add_argument("sequence_1", type=str)
parser.add_argument("sequence_2", type=str)

args = parser.parse_args()

# e.g. $ python main.py dynprog "ABCD" "[[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]]" "AAAAACCDDCCDDAAAAACC" "CCAAADDAAAACCAAADDCCAAAA"

if __name__ == "__main__":
    args.matrix = json.loads(args.matrix)
    if args.mode == "dynprog":
        result = local_alignment.dynprog(args.alphabet, args.matrix, args.sequence_1, args.sequence_2)
    elif args.mode == "dynproglin":
        result = local_alignment.dynproglin(args.alphabet, args.matrix, args.sequence_1, args.sequence_2)
    elif args.mode == "heuralign":
        result = local_alignment.heuralign(args.alphabet, args.matrix, args.sequence_1, args.sequence_2)
    else:
        print("Please specify execution mode mode: (dynprog / dynproglin / heuralign).")
        exit()
    print("Score: " + str(result[0]))
    print("Indices: " +str(result[1]) + str(result[2]))
