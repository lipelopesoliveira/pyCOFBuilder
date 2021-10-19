import os

COFs = os.listdir()

COFs.remove('original')
COFs.remove('rodando')

term = 0
fila = 0
erro = 0

for cof in COFs:
	out_file = os.path.join(os.getcwd(), cof, 'OPT_1', cof + '.cell_opt.out')
	try:
		tmp = open(out_file, 'r').readlines()
		terminou = False
		for i in tmp:
			if 'GEOMETRY OPTIMIZATION COMPLETED' in i:
				print(cof)
				terminou = True
				term += 1
				os.system(f'rm {cof}/OPT_1/*.wfn*')
				os.system(f'rm {cof}/OPT_1/*.restart')
		if terminou == False:
			print(cof, 'NAO TERMINOU!!!!!')
			erro += 1

	except:
		print('Nao achei o arquivo .out')
		fila += 1 
total = erro + term + fila

print(f'{term} finalizados. {term/total*100:.3f}%')
print(f'{fila} na fila. {fila/total*100:.3f}%')
print(f'{erro} deram erro {erro/total*100:.3f}%')
