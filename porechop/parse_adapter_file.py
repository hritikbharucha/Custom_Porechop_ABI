from .adapters import Adapter
import re
import sys

# Custom adapters have to be sumbmited as DNA.
DNA_rule = re.compile(r"^[ATCG]+$")

def get_adapters(infile):
    """ Parse a formatted text file containing custom adapter sequences
    Format example:
    line 1: Adapter name
    line 2: Start adapter sequence
    line 3: End adapter sequence
    --- repeat ---

    If you don't want to include the start or end sequence,
    just put an empty line.

    Each adapters will be named using the "adapter name" as a prefix, like:
    "adapter_name_Top"   : sequence
    "adapter_name_Bottom": sequence

    @param infile: path to the file containing the custom adapters
    @return adapters: list of Adapter object build from the input file.
    """
    adapters = []
    with open(infile) as f:
        loop = True
        while(loop):
            # try to get name and star/end sequences
            try:
                name = next(f).rstrip("\n")
                start = next(f).rstrip("\n")
                end = next(f).rstrip("\n")
            except StopIteration:
                loop = False
            else:
                ok = True
                # Check the sequences
                if(start):
                    start_test = DNA_rule.match(start)
                    ok = start_test is not None
                if(end):
                    end_test = DNA_rule.match(end)
                    ok = end_test is not None
                
                # If sequences are not in DNA format, quit.
                if(not ok):
                    print("/!\\\t Some sequences in the custom adapter file are not DNA.")
                    print("Start: ", start)
                    print("End: ", end)
                    exit(1)
                # Else, add everything to the adapter list.
                start_pair = (f"{name}_Top",start) if start else []
                end_pair = (f"{name}_Bottom",end) if end else []
                adapters.append(Adapter(name,
                                        start_sequence=start_pair,
                                        end_sequence=end_pair))
    return(adapters)