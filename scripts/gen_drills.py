# To run this script:
# - launch KiCAD's PCB editor
# - open Tools > Scripting Console
# - type "import gen_drills"

import os
import pcbnew

drills_file = 'outputs/drills.txt'
board_file = 'outputs/drills.kicad_pcb'

def create_via(board, net, x, y, drillsize, width):
    via = pcbnew.PCB_VIA(board)
    board.Add(via)
    via.SetNet(net)
    via.SetPosition(pcbnew.VECTOR2I_MM(x, y))
    via.SetDrill(int(drillsize * pcbnew.PCB_IU_PER_MM))
    via.SetWidth(int(width * pcbnew.PCB_IU_PER_MM))
    via.SetLayerPair(
        board.GetLayerID('F.Cu'),
        board.GetLayerID('B.Cu')
    )
    via.SetViaType(pcbnew.VIATYPE_THROUGH)

board = pcbnew.NewBoard(board_file)
ground = pcbnew.NETINFO_ITEM(board, "GND")
board.Add(ground)
for line in open(drills_file):
    x, y = map(float, line.split())
    create_via(board, ground, x, y, 0.5, 1.0)

pcbnew.SaveBoard(board_file, board)

