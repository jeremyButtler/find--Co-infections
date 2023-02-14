file ./cluster
b clusterGraph.c:143
b 225
b 251
b 260
b 289
r -f "testCodeFiles/longTest.tsv"

define pNId
    set $intCnt = 0
    output $intCnt
    set $nodeOn = $arg0
    set $intCnt = 1
    while $nodeOn
        output $intCnt
        print $nodeOn -> seqIdNode -> idCStr
        set $nodeOn = $nodeOn -> nextNode
        set $intCnt = $intCnt + 1
    end
end

define pEId
    set $intCnt = 0
    print $arg0 -> seqIdNode -> idCStr
    output $intCnt
    set $nodeOn = $arg0 -> edgesList
    set $intCnt = 1
    while $nodeOn
        output $intCnt
        print $nodeOn->seqIdNode->idCStr
        set $nodeOn = $nodeOn->nextNode
        set $intCnt = $intCnt + 1
    end
end


define pNStruct
    set $nodeOn = $arg0
    while $nodeOn
        print *$nodeOn
        set $nodeOn = $nodeOn -> nextNode
    end
end

define pEStruct
    print *$arg0
    set $nodeOn = $arg0 -> edgesList
    while $nodeOn
        print *$nodeOn
        set $nodeOn = $nodeOn -> nextNode
    end
end
