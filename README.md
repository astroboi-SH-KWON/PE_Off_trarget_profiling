# PE_Off_trarget_profiling


1. PAM 찾기
    PAM은 NGG / NAG / NGA
    PAM 부분도 mismatch가 있다면 개수 세기

2. PAM 기준 왼쪽 20개, 오른쪽 RTT-6) 서열을 target reference와 비교
3. mismatch 개수 파악 (2 <= mismatch <= 10)



나와야 할 정보 예시					
Type	Chr	Location	Off-target sequence	mismatch_cnt	strand
HEK3-01_guide mismatch	1	1453414	GGCCCAGACTGAGCACGTGATGGCAGAGGA		


::: 20201112 :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

1. Total mismatch count 범위를 '0<= mismatch count <= 8' 로 수정 (기존: 2~10)
- Mismatch 0인 것은 무조건 1개가 나오겠지만 같이 보면 좋을 것 같습니다.
2. Target 영역별 mismatch를 나타내고 싶습니다. Target 영역은 3개입니다.
- Guide (앞에서부터 1~20번째)
- PAM (21~23번째)
- RTT only (24번부터 끝까지)
3. 결과가 다 나온 다음에는 모든 리스트를 하나의 엑셀파일로 합쳐주시기 바랍니다.

나와야 할 정보 예시					
Type	Chr	Location	Off-target sequence	mismatch_tot_cnt	mismatch_cnt_in_guide	mismatch_cnt_in_PAM	mismatch_cnt_in_RTTonly	strand

