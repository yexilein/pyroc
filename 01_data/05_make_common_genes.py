
import sys
sys.path.append("..")
import data


def main():
    make_mouse_genes()
    
def make_mouse_genes():
    go_sets = data.load_go_mouse()
    coexp = data.mouse_coexp()
    ppi = data.mouse_ppi()
    common_genes = set(go_sets[1]).intersection(coexp[1]).intersection(ppi[1])
    with open("data/mouse_common_genes.txt", "w") as fp:
        fp.write("\n".join(common_genes))
        fp.write("\n")


if __name__ == "__main__":
    main()
