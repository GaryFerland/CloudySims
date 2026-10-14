//-------------------------------------------------------------------------------
//  VectorHash - a very fast hash function optimized using SIMD instructions
//
//  Copyright (c) 2018-2026 Peter A.M. van Hoof
//  All Rights Reserved
//
//  Distributed under the "zlib license". See the accompanying LICENSE file.
//-------------------------------------------------------------------------------

#include "TestMain.h"
#include "vectorhash.h"

namespace {

	// test empty buffer
	TEST(TestZeroLengthBufferScalar_32)
	{
		VectorHashSIMD(buffer, 0, 0xfd4c799d, cksum, IS_SCALAR, 32);
		CHECK( CheckHash(cksum, "7647d9bd") );
	}

	TEST(TestZeroLengthBufferScalar_64)
	{
		VectorHashSIMD(buffer, 0, 0xfd4c799d, cksum, IS_SCALAR, 64);
		CHECK( CheckHash(cksum, "73711a77d6031b6f") );
	}

	TEST(TestZeroLengthBufferScalar_128)
	{
		VectorHashSIMD(buffer, 0, 0xfd4c799d, cksum, IS_SCALAR, 128);
		CHECK( CheckHash(cksum, "fe82e7d9998e9819c7ac954ea0a0ea8e") );
	}

	TEST(TestZeroLengthBufferScalar_256)
	{
		VectorHashSIMD(buffer, 0, 0xfd4c799d, cksum, IS_SCALAR, 256);
		CHECK( CheckHash(cksum, "c87638f19f1778d7ee03dddc5c6677fd0701e532d311bfffddfc21c0ece8910a") );
	}

	TEST(TestZeroLengthBufferScalar_512)
	{
		VectorHashSIMD(buffer, 0, 0xfd4c799d, cksum, IS_SCALAR, 512);
		CHECK( CheckHash(cksum, "280f647d93f843b75627cc8c437084e2668a3cd9e077a2cf3186c67c990983c2"
		                        "8ff8180187d362a79cb02cb49ba350a5e185e0050fc273e836b1db30571cfa66") );
	}

	TEST(TestZeroLengthBufferScalar_1024)
	{
		VectorHashSIMD(buffer, 0, 0xfd4c799d, cksum, IS_SCALAR, 1024);
		CHECK( CheckHash(cksum, "c38bd292e55400931e47f4397b5266b3e143ffd918fb944c33ad908458a08c62"
		                        "65e9c28007bf68301ab6f68c82ddcc2fa993ad69556b9ab4d588594ac4d02b5d"
		                        "dd083a4cf1418b6522f06bff90a67bf93faaa719e1dbdbda3592deddbd1a8378"
		                        "f0993653e73c292f338729c1efe01b352e90ed9a0ae01e613ad72736db619b2d") );
	}

	// test half a blocksize
	TEST(Test128ByteBufferScalar_32)
	{
		CHECK( ReadBuffer("test0128", 128, buffer) );
		VectorHashSIMD(buffer, 128, 0xfd4c799d, cksum, IS_SCALAR, 32);
		CHECK( CheckHash(cksum, "50c5164c") );
	}

	TEST(Test128ByteBufferScalar_64)
	{
		CHECK( ReadBuffer("test0128", 128, buffer) );
		VectorHashSIMD(buffer, 128, 0xfd4c799d, cksum, IS_SCALAR, 64);
		CHECK( CheckHash(cksum, "a9b41194bceedd04") );
	}

	TEST(Test128ByteBufferScalar_96)
	{
		CHECK( ReadBuffer("test0128", 128, buffer) );
		VectorHashSIMD(buffer, 128, 0xfd4c799d, cksum, IS_SCALAR, 96);
		CHECK( CheckHash(cksum, "b01e0af4bceedd04b71200b9") );
	}

	TEST(Test128ByteBufferScalar_128)
	{
		CHECK( ReadBuffer("test0128", 128, buffer) );
		VectorHashSIMD(buffer, 128, 0xfd4c799d, cksum, IS_SCALAR, 128);
		CHECK( CheckHash(cksum, "970566856108ce7bce658104386829fa") );
	}

	TEST(Test128ByteBufferScalar_160)
	{
		CHECK( ReadBuffer("test0128", 128, buffer) );
		VectorHashSIMD(buffer, 128, 0xfd4c799d, cksum, IS_SCALAR, 160);
		CHECK( CheckHash(cksum, "efa539560f142d7e5d1a3dc9b81f7bdf2c2b94f8") );
	}

	TEST(Test128ByteBufferScalar_192)
	{
		CHECK( ReadBuffer("test0128", 128, buffer) );
		VectorHashSIMD(buffer, 128, 0xfd4c799d, cksum, IS_SCALAR, 192);
		CHECK( CheckHash(cksum, "6a109d9852be15a9053bb48121de06bf95006085df2f9c55") );
	}

	TEST(Test128ByteBufferScalar_224)
	{
		CHECK( ReadBuffer("test0128", 128, buffer) );
		VectorHashSIMD(buffer, 128, 0xfd4c799d, cksum, IS_SCALAR, 224);
		CHECK( CheckHash(cksum, "ab99379bd8e63165917952a677dc830b673d34d9eecd69acbc2fcec7") );
	}

	TEST(Test128ByteBufferScalar_256)
	{
		CHECK( ReadBuffer("test0128", 128, buffer) );
		VectorHashSIMD(buffer, 128, 0xfd4c799d, cksum, IS_SCALAR, 256);
		CHECK( CheckHash(cksum, "630b825fec4275e28b02590c53db6d920e3647133d5300a1e2887861867dfcf0") );
	}

	TEST(Test256ByteBufferScalar_384)
	{
		CHECK( ReadBuffer("test0256", 256, buffer) );
		VectorHashSIMD(buffer, 256, 0xfd4c799d, cksum, IS_SCALAR, 384);
		CHECK( CheckHash(cksum, "89a883ef58dd0a9fe9ea10f66f950a9e57f9a12faaf96767a6d26a8b4d6b6bfd"
		                        "c3bd176e0ed69ebfe286c38acfad26e5") );
	}

	TEST(Test256ByteBufferScalar_512)
	{
		CHECK( ReadBuffer("test0256", 256, buffer) );
		VectorHashSIMD(buffer, 256, 0xfd4c799d, cksum, IS_SCALAR, 512);
		CHECK( CheckHash(cksum, "f0ef0c9efd9c02c38be2f2687ac20acee763fe6063829747aace159ddcdc7470"
		                        "a7e4cdbd00c7ed0fbc30c8810fd35b73f7b65e5242369a63d3ea6c2e7c0ae56e") );
	}

	TEST(Test512ByteBufferScalar_768)
	{
		CHECK( ReadBuffer("test0512", 512, buffer) );
		VectorHashSIMD(buffer, 512, 0xfd4c799d, cksum, IS_SCALAR, 768);
		CHECK( CheckHash(cksum, "bf781a572cd0fb74d6c22c92f501b7c2109541f079929ff08b6ace8956de1280"
		                        "4e7afa6bbfdcd9eb727e8766780ccb35b8f8392c59329aa8d93d5ab5149e0495"
		                        "73f97958e25805a5ea56b5a447643d4c707ef066577f3cab105bc3e9c6a96e2b") );
	}

	TEST(Test512ByteBufferScalar_1024)
	{
		CHECK( ReadBuffer("test0512", 512, buffer) );
		VectorHashSIMD(buffer, 512, 0xfd4c799d, cksum, IS_SCALAR, 1024);
		CHECK( CheckHash(cksum, "4f3d9728dfadf99677de56e642fc5990e3133fe2008cf431c93931011d66f373"
		                        "4d9ed038d6d581f5e8bf073cec644d564455bbae6c70c9a5d1bd3631835478b2"
		                        "bec773a1c1daecd260ca92cdaade67869c03e4fc5cd1f14d50addab0d34aa534"
		                        "5e1b103946534375900f08abd8195b9bf22592eb44bbc1fc37eb1d6e5039014e") );
	}

	// test one blocksize
	TEST(Test256ByteBufferScalar_32)
	{
		CHECK( ReadBuffer("test0256", 256, buffer) );
		VectorHashSIMD(buffer, 256, 0xfd4c799d, cksum, IS_SCALAR, 32);
		CHECK( CheckHash(cksum, "0eabedeb") );
	}

	TEST(Test256ByteBufferScalar_64)
	{
		CHECK( ReadBuffer("test0256", 256, buffer) );
		VectorHashSIMD(buffer, 256, 0xfd4c799d, cksum, IS_SCALAR, 64);
		CHECK( CheckHash(cksum, "81534826f50039f6") );
	}

	TEST(Test256ByteBufferScalar_96)
	{
		CHECK( ReadBuffer("test0256", 256, buffer) );
		VectorHashSIMD(buffer, 256, 0xfd4c799d, cksum, IS_SCALAR, 96);
		CHECK( CheckHash(cksum, "99a9d7b4f50039f60c5afc8d") );
	}

	TEST(Test256ByteBufferScalar_128)
	{
		CHECK( ReadBuffer("test0256", 256, buffer) );
		VectorHashSIMD(buffer, 256, 0xfd4c799d, cksum, IS_SCALAR, 128);
		CHECK( CheckHash(cksum, "b384fcaef98a474cd776074b9d78bca9") );
	}

	TEST(Test256ByteBufferScalar_160)
	{
		CHECK( ReadBuffer("test0256", 256, buffer) );
		VectorHashSIMD(buffer, 256, 0xfd4c799d, cksum, IS_SCALAR, 160);
		CHECK( CheckHash(cksum, "adf261ccc63e8066e51e6c096d32b39fcddbf78e") );
	}

	TEST(Test256ByteBufferScalar_192)
	{
		CHECK( ReadBuffer("test0256", 256, buffer) );
		VectorHashSIMD(buffer, 256, 0xfd4c799d, cksum, IS_SCALAR, 192);
		CHECK( CheckHash(cksum, "7d3a49576f3a7d8799a7bdfeef254e39a8712c1bcbc5ad54") );
	}

	TEST(Test256ByteBufferScalar_224)
	{
		CHECK( ReadBuffer("test0256", 256, buffer) );
		VectorHashSIMD(buffer, 256, 0xfd4c799d, cksum, IS_SCALAR, 224);
		CHECK( CheckHash(cksum, "32381fd2605a5710f0a2c49f1194dfc036e6d105e82130c1e91f57a9") );
	}

	TEST(Test256ByteBufferScalar_256)
	{
		CHECK( ReadBuffer("test0256", 256, buffer) );
		VectorHashSIMD(buffer, 256, 0xfd4c799d, cksum, IS_SCALAR, 256);
		CHECK( CheckHash(cksum, "96750ccd5f66ee6fcf17f512abecd3339e487a437976ad5947589e070d8e8d9e") );
	}

	TEST(Test512ByteBufferScalar_384)
	{
		CHECK( ReadBuffer("test0512", 512, buffer) );
		VectorHashSIMD(buffer, 512, 0xfd4c799d, cksum, IS_SCALAR, 384);
		CHECK( CheckHash(cksum, "778622bc016f4e7d40d4a2ad75a1b244938c9d03948560d067bf37e935c11ccd"
		                        "0ade4907afbf792a384307dc67e9d867") );
	}

	TEST(Test512ByteBufferScalar_512)
	{
		CHECK( ReadBuffer("test0512", 512, buffer) );
		VectorHashSIMD(buffer, 512, 0xfd4c799d, cksum, IS_SCALAR, 512);
		CHECK( CheckHash(cksum, "e6ad8e7d595ad7e025cce2f56aea6fdc01c75ac38d68234a6f4b0ce5fb335dcb"
		                        "7790464c11cb5f00e5c15559911784f11f12ac10a75ccd3208f0f99e4a35ac4b") );
	}

	TEST(Test1024ByteBufferScalar_768)
	{
		CHECK( ReadBuffer("test1024", 1024, buffer) );
		VectorHashSIMD(buffer, 1024, 0xfd4c799d, cksum, IS_SCALAR, 768);
		CHECK( CheckHash(cksum, "5d89a5a845138e4e381349667a88d9a85eab9e32872214602dceff17d9558512"
		                        "6d25845b7c37c845d81242ba7f05fb9c2bbce954bc17a3c730f374df4d9948b8"
		                        "fcee4c9d13585b3e7626a59bf50b0339f27095296f1c0899238c5f3f41358d7f") );
	}

	TEST(Test1024ByteBufferScalar_1024)
	{
		CHECK( ReadBuffer("test1024", 1024, buffer) );
		VectorHashSIMD(buffer, 1024, 0xfd4c799d, cksum, IS_SCALAR, 1024);
		CHECK( CheckHash(cksum, "7cb508d7ccd3d6d7d025af8bfa7e6bafdad15d789ff09b0ee3e4e16de58ec936"
		                        "01ece64099f0273ab5f847c6907fb0a65ae1571a6a4a17ca98376173a71c6635"
		                        "8f8886bb85cb56185a6a09b315cd1a9494b68a913bb15d4a969a2312f801af24"
		                        "c6f293256f2feb733d957a3f95a195e2cb91b6202e321258a653a94cad0c8753") );
	}

	// test 3/2 blocksizes
	TEST(Test384ByteBufferScalar_32)
	{
		CHECK( ReadBuffer("test0384", 384, buffer) );
		VectorHashSIMD(buffer, 384, 0xfd4c799d, cksum, IS_SCALAR, 32);
		CHECK( CheckHash(cksum, "7273a9bf") );
	}

	TEST(Test384ByteBufferScalar_64)
	{
		CHECK( ReadBuffer("test0384", 384, buffer) );
		VectorHashSIMD(buffer, 384, 0xfd4c799d, cksum, IS_SCALAR, 64);
		CHECK( CheckHash(cksum, "e65f645584bcb59a") );
	}

	TEST(Test384ByteBufferScalar_128)
	{
		CHECK( ReadBuffer("test0384", 384, buffer) );
		VectorHashSIMD(buffer, 384, 0xfd4c799d, cksum, IS_SCALAR, 128);
		CHECK( CheckHash(cksum, "00d2ad3ff87e9969156227bcedce13ea") );
	}

	TEST(Test384ByteBufferScalar_256)
	{
		CHECK( ReadBuffer("test0384", 384, buffer) );
		VectorHashSIMD(buffer, 384, 0xfd4c799d, cksum, IS_SCALAR, 256);
		CHECK( CheckHash(cksum, "6838a170786ba5b21c0bce56feda9028749a294162501d660ad4ecdeee9c8245") );
	}

	TEST(Test768ByteBufferScalar_512)
	{
		CHECK( ReadBuffer("test0768", 768, buffer) );
		VectorHashSIMD(buffer, 768, 0xfd4c799d, cksum, IS_SCALAR, 512);
		CHECK( CheckHash(cksum, "57c8793d85812151b68a835a9ae6654b08d1ccc59a70118f4d6a6bd46be0e3f7"
		                        "ae3c9f83fab4d27dd9254a80ccde663d43a9ed48000f82f5b8fefcbdf0251957") );
	}

	TEST(Test1536ByteBufferScalar_1024)
	{
		CHECK( ReadBuffer("test1536", 1536, buffer) );
		VectorHashSIMD(buffer, 1536, 0xfd4c799d, cksum, IS_SCALAR, 1024);
		CHECK( CheckHash(cksum, "9377ee81b5ee4163f6213f395df47551a46ae82ca78b8bcd5c0e49c0606b7136"
		                        "be2dfb4949aa69e59d2ab4809cdc3b4676d133ca8a8ef884c88e93ce6add5875"
		                        "1fe5a21aeb0375d98299cff1634fa10bb4032a4912fe6ea18487c6e9eaab7cb0"
		                        "e0c9edbb5276ae17bf089013eb027c459db57c491d9bba404e365787eff9dafe") );
	}

	// test a large buffer
	TEST(TestLargeBufferScalar_32)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0xfd4c799d, cksum, IS_SCALAR, 32);
		CHECK( CheckHash(cksum, "81dddc5e") );
	}

	TEST(TestLargeBufferScalar_64)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0xfd4c799d, cksum, IS_SCALAR, 64);
		CHECK( CheckHash(cksum, "021fa7d3f643f834") );
	}

	TEST(TestLargeBufferScalar_128)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0xfd4c799d, cksum, IS_SCALAR, 128);
		CHECK( CheckHash(cksum, "5c3c9fb8481be32ea676886ab251f4fc") );
	}

	TEST(TestLargeBufferScalar_256)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0xfd4c799d, cksum, IS_SCALAR, 256);
		CHECK( CheckHash(cksum, "888a7ef5d0feda6a571ca947ece56c8833019cc8dca62dad01e6f9e60dcc29d3") );
	}

	TEST(TestLargeBufferScalar_512)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0xfd4c799d, cksum, IS_SCALAR, 512);
		CHECK( CheckHash(cksum, "147c905b96a3e78f458bb56cbab1f17c843f70d485b24bfa51d1823602481eea"
		                        "7890edaf1575a8ef04a57c92fc86353ff5c48c3c8198762f6b409d9aa52bff52") );
	}

	TEST(TestLargeBufferScalar_1024)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0xfd4c799d, cksum, IS_SCALAR, 1024);
		CHECK( CheckHash(cksum, "3dfe155d7c8fab0ccf04baf748bf2bd9a0dd00924eb7d3ba94af78ba1976232f"
		                        "a4ec9c46080893e4321601b07f0e7aeac6a4ab010e505d443889fc978b23b0f9"
		                        "44ca22d4e4f6079b94b8985dd1adc1b13f4d2a2c08dd53e22f3bb84eb0fbab75"
		                        "06eff60a3089394a57f98b532a89cd66c7d6d2da971025c0cad83c09a3aa3d21") );
	}

	// test alternative seed
	TEST(TestLargeBufferAltSeedScalar_32)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0x6ec74615, cksum, IS_SCALAR, 32);
		CHECK( CheckHash(cksum, "2717e727") );
	}

	TEST(TestLargeBufferAltSeedScalar_64)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0x6ec74615, cksum, IS_SCALAR, 64);
		CHECK( CheckHash(cksum, "e00c65eadc62f6d4") );
	}

	TEST(TestLargeBufferAltSeedScalar_128)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0x6ec74615, cksum, IS_SCALAR, 128);
		CHECK( CheckHash(cksum, "b26873b350a8095d2ded60d9cf2d1a37") );
	}

	TEST(TestLargeBufferAltSeedScalar_256)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0x6ec74615, cksum, IS_SCALAR, 256);
		CHECK( CheckHash(cksum, "79c87bd2c4d5368c277c9fc2b1b871f0ec6077b3a129a62bdcd4bee3ba44cc17") );
	}

	TEST(TestLargeBufferAltSeedScalar_512)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0x6ec74615, cksum, IS_SCALAR, 512);
		CHECK( CheckHash(cksum, "4661b2e71c540ad70f9234124d20482dcf891cbd727933cb9a574448d66407aa"
		                        "e4b3f7b34c200e84efd84a82ab6f26bf8796dd57ea89503bb62369dcde5cd921") );
	}

	TEST(TestLargeBufferAltSeedScalar_1024)
	{
		CHECK( ReadBuffer("test9999", 1048576, buffer) );
		VectorHashSIMD(buffer, 1048576, 0x6ec74615, cksum, IS_SCALAR, 1024);
		CHECK( CheckHash(cksum, "6874b6e17181cade7bef55f316ee4df8bfb48131baca25484f97bf5234650fc3"
		                        "4873d6ae572d07fa01df034e45a62db004f32f7871427c0db41dd758ee32f0c3"
		                        "64bbd5c252b0d96e658c2c183202d8528530a57fc21dc074274f2e24e0abdb2a"
		                        "8126b974494c396114937f68312a076c898554e460dda37fc5e103ef5ee79f1e") );
	}

}
